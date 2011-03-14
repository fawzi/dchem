/// calculator that gets contexts from external calculators that announce themselves to this calculator
module dchem.calculator.CalculatorServer;
import dchem.calculator.CalculatorModels;
import dchem.calculator.Calculator;
import dchem.sys.ParticleSys;
import dchem.input.RootInput;

// at the moment very simple: eagerly gets all contexts
// should probably wrap the contexts and be smarter and more aware of their use (supporting history for example)

interface CalculatorServerI{
    void addContextWithUrl(char[]url);
    //void addClientWithUrl(char[]url);
}

interface CalculatorClientI{
    void addAllContextsTo(char[]url);
}

class CalculatorServerGen:Method{
    InputField startClient;
    char[] logFile;
    
    mixin(serializeSome("dchem.CalculatorServer",`calculator that sets up a server to which real calculators can connect`,`
    startClient: client started when the first request comes in
    logFile: place where to log the url that can be used to access this server`));
    mixin printOut!();
    mixin myFieldMixin!();
    bool verify(CharSink s){
        bool res=true;
        if (startClient !is null){
            if ((cast(CalculatorClientI)startClient.contentObj)is null){
                res=false;
                dumper(s)("startClient for field ")(myFieldName)(" if given should be a CalculatorClient\n");
            }
        }
        return res;
    }
    
    Deque!(CalculationContext) readyContexts;
    WaitCondition waitContexts;
    bool didStartClient;
    //HashTable!(char[],CalculationClient) availableClients;
    //Deque!(CalculationClient) readyClients;
    
    /// adds the context with the given url
    void addContextWithUrl(char[]url){
        auto ctx=ProtocolHandler.proxyForUrlT!(CalculationContext)(url);
        readyContexts.pushBack(ctx);
        waitContexts.checkCondition();
    }
    /// gets a calculator to perform calculations with this method
    CalculationContext getCalculator(bool wait,ubyte[]history){
        if (!didStartClient){
            synchronized(this){
                if (!didStartClient){
                    didStartClient=true;
                    if (startClient!is null && (cast(CalculatorClientI)startClient.contentObj)!is null){
                        auto client=cast(CalculatorClientI)startClient.contentObj;
                        client.addAllContextsTo(externalUrl);
                    }
                }
            }
        }
        CalculationContext res;
        while(true){
            if (readyContexts.popBack(res)) return res;
            if (!wait) return null;
            waitContexts.wait();
        }
    }
    /// drops the history associated with the given key
    void dropHistory(ubyte[]history){}
    /// clears all history
    void clearHistory(){}
    mixin(rpcMixin("dchem.CalculatorServer","CalculatorServerI",`
    addContextWithUrl`));
    DefaultVendor vendor;
    this(){
        vendor=new DefaultVendor(this);
        ProtocolHandler.defaultProtocol.publisher.publishObject(vendor,"calculatorClient_"~myFieldName);
    }
}

/// generates the 
class CalculatorClientGen:Sampler,CalculatorClientI{
    char[] serverUrl;
    InputField method;
    
    mixin myFieldMixin!();
    bool verify(CharSink s){
        bool res=true;
        if (method is null || (cast(Method)method.contentObj)is null){
            res=false;
            dumper(s)("method for field ")(myFieldName)(" should be a valid method\n");
        }
        return res;
    }
    
    void addAllContextsTo(char[]url){
        serverUrl=url;
        auto server=ProtocolHandler.proxyForUrlT!(CalculatorServerI)(url);
        while (true){
            ctx=getCalculator(true,null);
            if (ctx is null) return null;
            server.addContextWithUrl(ctx.externalUrl);
        }
    }
    
    void run(){
        addAllContextsTo(serverUrl);
    }
    mixin(serializeSome("dchem.CalculatorClient",`defines a client that will connect to a CalculatorServer`,`
    serverUrl: the url to which this client should connect
    method: the method to use for the context creation`));
    mixin printOut!();
    mixin(rpcMixin("dchem.CalculatorClient", "CalculatorClientI","addAllContextsTo:oneway"));
    
    DefaultVendor vendor;
    this(){
        vendor=new DefaultVendor(this);
        ProtocolHandler.defaultProtocol.publisher.publishObject(vendor,"calculatorClient_"~myFieldName);
    }
    char[]exportedUrl(){
        return mainVendor.proxyObjUrl();
    }
}