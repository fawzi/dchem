module dchem.calculator.Cp2kCalculator;
import dchem.calculator.Calculator;
import dchem.calculator.FileCalculator;
import dchem.sys.ParticleSys;
import dchem.input.RootInput;
import blip.parallel.mpi.MpiModels;
import Search=tango.text.Search;
import tango.sys.Process;
import dchem.calculator.ProcContext;
import blip.io.Socket;
import blip.io.BasicIO;
import blip.io.BufferIn;
import blip.io.BasicStreams;
import blip.parallel.smp.Smp;
import blip.io.EventWatcher;
import blip.serialization.Serialization;
import blip.container.Deque;
import blip.sync.Atomic;
import blip.io.Console;
import blip.math.random.Random;
import blip.container.GrowableArray;
import dchem.util.ExecCmd;
import blip.io.SimpleBinaryProtocol;
import blip.bindings.blas.Types;
import blip.container.BulkArray;
import blip.math.Math;
import dchem.Common;

/// represent a connection with an host
class Cp2kConnection{
    enum Status{
        Setup,
        Running,
        Stopping,
        Stopped,
    }
    Cp2kServer server;
    BasicSocket sock;
    BufferedBinStream outStream;
    BufferIn!(void) readIn;
    Reader!(char) charReader;
    SequentialTask commTask; /// sequential task that should be used to communicate on this connection
    size_t localUsers=1;
    TargetHost targetHost;
    Status status=Status.Setup;
    CharSink log;
    LoopHandlerI loop;
    ev_tstamp created; // use this to kill connections that have to wait too much??
    ev_tstamp lastUse;
    char[] connectionId;
    
    override equals_t opEquals(Object o){
        return this is o;
    }
    override int opCmp(Object o){
        size_t a=cast(size_t)cast(void*)this;
        size_t b=cast(size_t)cast(void*)o;
        return ((a<b)?-1:((a==b)?0:1));
    }
    final void writeExact(void[] src){
        this.sock.writeExactTout(src,loop);
    }
    final size_t rawReadInto(void[] dest){
        return this.sock.rawReadIntoTout(dest,loop);
    }
    final void rawReadExact(void[]dest){
        readExact(&this.rawReadInto,dest);
    }
    this(Cp2kServer cp2kServ,TargetHost targetHost,BasicSocket sock){
        this.server=cp2kServ;
        this.targetHost=targetHost;
        this.sock=sock;
        this.loop=cp2kServ.loop;
        this.created=loop.now();
        version(Cp2kNoCache){}
        else {
            this.sock.noDelay(true); // use no delay to reduce the latency
        }
        //this.sock.keepalive(true);
        commTask=new SequentialTask("cp2kCommTask",defaultTask,true);
        // should limit buffer to 1280 or 1500 or multiples of them? (jumbo frames)
        outStream=new BufferedBinStream("cp2kServer",&this.writeExact,3000,&this.sock.flush,&this.sock.close);
        readIn=new BufferIn!(void)(&outStream.desc,&this.rawReadInto);
        localUsers=1;
        log=cp2kServ.log;
        version(TrackCp2kSock){
            sinkTogether(log,delegate void(CharSink s){
                dumper(s)("created new cp2k connection to ")(targetHost)(" on socket ")(this.sock)("\n");
            });
        }
        char[128] buf;
        connectionId=recv(buf,false).dup;
        version(TrackCp2kSock){
            sinkTogether(log,delegate void(CharSink s){
                dumper(s)("cp2k connection to ")(targetHost)(" on socket ")(this.sock)(" has id '")(connectionId)("'\n");
            });
        }
    }
    
    this(Cp2kServer cp2kServ,TargetHost targetHost){
        this(cp2kServ,targetHost,BasicSocket(targetHost));
    }
    
    bool tryAddLocalUser(){
        synchronized(this){
            if (localUsers==0) return false;
            ++localUsers;
            return true;
        }
    }

    // stuff to keep track of the users (and allow in the future to close a unused connection)
    void addLocalUser(){
        synchronized(this){
            if (atomicAdd(localUsers,1)==0){
                throw new Exception("localUsers was 0 in addLocalUser",__FILE__,__LINE__);
            }
        }
    }
    void rmLocalUser(){
        synchronized(this){
            auto oldL=atomicAdd(localUsers,-1);
            if (oldL==0){
                throw new Exception("localUsers was 0 in rmLocalUser",__FILE__,__LINE__);
            }
            if (oldL==1){
                lastUse=loop.now;
            }
        }
    }
    /// sends a sbp type (char[],int[],double[])
    void send(T)(T arr){
        sbpSend(&writeExact,arr);
    }
    /// receives a sbp type (char[],int[],double[])
    DynamicArrayType!(T) recv(T)(T buf,bool strict=true){
        return sbpRead(&rawReadExact,buf,strict);
    }
    // closes the connection
    void closeConnection(){
        synchronized(this){
            if (status==Status.Running){
                status=Status.Stopping;
            } else {
                status=Status.Stopped;
            }
        }
        outStream.close();
        readIn.shutdownInput();
    }
}

class Cp2kServer:InputElement{
    ushort portMin=49152;
    ushort portMax=65535;
    char[] port="51000";
    char[] setupCmd;
    double maxLife=525600.0;
    bool strict=false;
    
    mixin(serializeSome("dchem.Cp2kServer",`
    portMin: the minimum port number to use as fallback (49152)
    portMax: the maximum port number to use as fallback (65535)
    port: the port/protocol that should be used if possible (51000)
    setupCmd: the command executed on activate (might start the clients, in it "[port]" is replaced with the port of the server)
    maxLife: the maximum life of a connection in minutes (525600 = one year)
    strict:if the port selection should be changed (false)`));
    mixin myFieldMixin!();
    mixin printOut!();
    
    SocketServer server;
    CharSink log;
    Deque!(Cp2kConnection) connections;
    WaitCondition waitConnection;
    LinearComm pEnv;
    enum Status{
        Setup,Starting,Running,Stopping
    }
    Status status;
    LoopHandlerI loop;
    
    void giveBackConnection(Cp2kConnection c){
        if (loop.now()>c.created+maxLife*60){
            c.closeConnection();
        } else {
            connections.pushBack(c);
            waitConnection.checkCondition();
        }
    }
    void setup(LinearComm pEnv,CharSink log){
        bool start=false;
        synchronized(this){
            if (status==Status.Setup){
                status=Status.Starting;
                pEnv=pEnv;
                start=true;
            } else {
                if (pEnv !is this.pEnv) throw new Exception("different pEnv",__FILE__,__LINE__);
                if (status>Status.Running) throw new Exception("setup called on stopped Cp2kServer",__FILE__,__LINE__);
                return;
            }
        }
        if (pEnv.myRank!=0) assert(0,"unimplemented");
        if (start){
            bool isBound=false;
            this.log=sout.call;
            server=new SocketServer(port,&this.handleConnection,log);
            Exception bindE;
            for (int i=0;i<100;++i){
                try{
                    server.start();
                    isBound=true;
                } catch(BIONoBindException e){
                    bindE=e;
                }
                if (isBound || strict) break;
                if (portMin>=portMax) throw new Exception("could not bind to requested port, and no fallback for Cp2kServer from filed "~myFieldName,__FILE__,__LINE__);
                auto newP=rand.uniformR2(portMin,portMax);
                char[25] buf;
                auto arr=lGrowableArray(buf);
                writeOut(&arr.appendArr,newP);
                port=arr.takeData();
                server.serviceName=port;
            }
            if (!isBound) {
                server=null;
                synchronized(this){
                    status=Status.Stopping;
                }
                throw new BIONoBindException("could not bind Cp2kServer from field "~myFieldName,__FILE__,__LINE__,bindE);
            }
            synchronized(this){
                status=Status.Running;
            }
            log("started Cp2kServer for field "~myFieldName~" on port "~port~"\n");
            if (setupCmd.length>0){
                auto repl=Search.find("[port]");
                auto realCmd=repl.replace(setupCmd,port);
                
                auto opInProgress=new Process(realCmd); // set workdir???
                int status;
                sinkTogether(log,delegate void(CharSink s){
                    dumper(s)("executing command:")(realCmd)("\n");
                });
                log(getOutput(opInProgress,status));
                if (status!=0){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("command "~realCmd~" failed with status ")(status);
                    }),__FILE__,__LINE__);
                }
            }
        }
    }
    
    Cp2kConnection getConnection(bool wait){
        while (true){
            Cp2kConnection newC;
            if (connections.popFront(newC)){
                if (status==Status.Stopping || newC.created+60*maxLife<loop.now()){
                    newC.closeConnection();
                    continue; // close all connections on stop...
                }
                return newC;
            }
            if (!wait || status==Status.Stopping) return null;
            waitConnection.wait();
        }
    }
    
    void handleConnection(ref SocketServer.Handler h){
        TargetHost th=h.otherHost();
        auto newC=new Cp2kConnection(this,th,h.sock);
        version(TrackRpc){
            sinkTogether(log,delegate void(CharSink s){
                dumper(s)(taskAtt.val)(" got a cp2k connection to ")(th)("\n");
            });
        }
        connections.pushBack(newC);
        waitConnection.checkCondition();
    }
    
    bool hasConnections(){
        return connections.length!=0 || status==Status.Stopping;
    }
    
    this(){
        waitConnection=new WaitCondition(&this.hasConnections);
        loop=noToutWatcher; // retain???
        log=sout.call;
    }
    
    bool verify(CharSink log){
        bool res=true;
        if (maxLife<=0){
            res=false;
            sinkTogether(log,delegate void(CharSink s){
                dumper(s)("maxLife in field ")(myFieldName)(" must be larger than 0\n");
            });
        }
        return res;
    }
}

class Cp2kMethod:TemplateExecuter{
    char[] initialFileToLoad;
    InputField cp2kServer;
    
    mixin(serializeSome("dchem.Cp2kMethod",`
    initialFileToLoad: the file to load when the setup is complete
    cp2kServer: the cp2k server to use to get connections to cp2k instances`));
    
    bool verify(CharSink log){
        return true;
    }
    Cp2kServer serv(){
        auto res=cast(Cp2kServer)cp2kServer.contentObj;
        assert(res!is null);
        return res;
    }
    void setup(LinearComm pWorld,CharSink log){
        serv.setup(pWorld,log);
        subs["[port]"]=serv.port;
        super.setup(pWorld,log);
    }
    bool verify(CharSink log){
        bool res=true;
        if ((cast(Cp2kServer)cp2kServer.contentObj)is null){
            res=false;
            sinkTogether(log,delegate void(CharSink s){
                dumper(s)("cp2kServer in field ")(myFieldName)(" must be set and of type Cp2kServer\n");
            });
        }
        return res;
    }
    CalculationContext getCalculator(bool wait,ubyte[]history){
        auto conn=serv.getConnection(wait);
        if (conn!is null){
            auto ctx=new Cp2kContext(conn,this,collectAppender(delegate void(CharSink s){
                dumper(s)("cp2k")(ProcContext.instance.id)("-")(ProcContext.instance.localId.next());
            }));
            return ctx;
        }
        return null;
    }
    
}

class Cp2kContext: ExecuterContext{
    Cp2kMethod cp2kM;
    Cp2kConnection connection;
    int[1] myCtxId;
    
    void checkReady(){
        char[128] buf;
        auto res=connection.recv(buf,false);
        if (res!="* READY") throw new Exception("unexpected reply:"~res,__FILE__,__LINE__);
    }

    this(Cp2kConnection connection,Cp2kMethod input,char[] contextId,TemplateHandler th=null){
        super(input,contextId,th);
        this.connection=connection;
        connection.send("LOAD "~cp2kM.initialFileToLoad);
        int [1] cId;
        myCtxId[]=connection.recv(myCtxId);
        checkReady();
        connection.send("GET_NATOM");
        connection.send(myCtxId);
        int[1] nat;
        connection.recv(nat);
        if (nat[0]!=pSysReal.dynVars.x.pos.dataLength) throw new Exception(collectAppender(delegate void(CharSink s){
            dumper(s)("unexpected number of atoms ")(nat[0])("vs")(pSysReal.dynVars.x.pos.dataLength);
        }));
    }
    
    void updateEF(bool updateE=true,bool updateF=true){
        if (updateE || updateF){
            checkReady();
            connection.send("SET_POS");
            connection.send(myCtxId);
            auto pos=pSysReal.dynVars.x.pos;
            auto externalOrder=pSysReal.sysStruct.externalOrder;
            auto toSend=BulkArray!(f_double)(
                
externalOrder.nLocalParticles(pSysReal.sysStruct.levels[0])*3+1);
            size_t ii=0;
            foreach (idx;externalOrder.gSortedLocalPIndex.sLoop){
                auto p=pos[idx,0];
                toSend[++ii]=p.x;
                toSend[++ii]=p.y;
                toSend[++ii]=p.z;
            }
            connection.send(toSend.data[1..$]);
            checkReady();
            if (updateE && updateF){
                connection.send("UPDATE_EF");
            } else if (updateE){
                connection.send("UPDATE_E");
            } else {
                connection.send("UPDATE_F");
            }
            connection.send(myCtxId);
            if (updateE){
                checkReady();
                connection.send("GET_E");
                double[2] e;
                connection.recv(e);
                pSysReal.dynVars.potentialEnergy=e[0];
                pSysReal.dynVars.potentialEnergyError=e[1];
            }
            if (updateF){
                checkReady();
                connection.send("GET_F");
                pSysReal.checkMddx();
                connection.recv(toSend.data);
                size_t ij=0;
                pSysReal.dynVars.mddxError=toSend[0];
                foreach (idx;externalOrder.gSortedLocalPIndex.sLoop){
                    Vector!(Real,3) p;
                    p.x=toSend[++ij];
                    p.y=toSend[++ij];
                    p.z=toSend[++ij];
                    pos[idx,0]=p;
                }
            }
            toSend.guard.release();
        }
        maxChange=0;
        changeLevelSet=ChangeLevel.SmoothPosChange;
    }
}
