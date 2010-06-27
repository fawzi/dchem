/// requested properties for the current evaluation
module dchem.sys.Requests;
import blip.core.Variant;
import blip.serialization.Serialization;

struct EnergyRequests{
    bool energy;
    bool force;
    int high_e_deriv;
    mixin(serializeSome("dchem.sys.Requests.EnergyRequests",
        `energy: the energy of the system
        force: derivatives with respect to position
        high_e_deriv: higher derivatives, meaningful only if larger than 1`));
}

class Requests{
    EnergyRequests energy;
    Variant[char[]] others;
}