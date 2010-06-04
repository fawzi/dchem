module dchem.sampler.MinEExplorer;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import tango.math.random.Random;
import blip.rtest.RTest;
import Atomic=blip.sync.Atomic;
import blip.t.math.Math:max;
import blip.container.BatchedGrowableArray;
import tango.util.container.HashSet;
import dchem.input.RootInput;
import blip.container.Deque;


alias ulong EKey; /// key of a MinEExplorer instance, 0 is an invalid key



alias HashSet Set;
