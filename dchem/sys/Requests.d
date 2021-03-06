/// requested properties for the current evaluation
module dchem.sys.Requests;
import blip.core.Boxer;
import blip.serialization.Serialization;

struct EnergyRequests{
    bool energy;
    bool force;
    int high_e_deriv;
    mixin(serializeSome("dchem.sys.Requests.EnergyRequests","pippo delete this file",
        `energy: the energy of the system
        force: derivatives with respect to position
        high_e_deriv: higher derivatives, meaningful only if larger than 1`));
}

class Requests{
    EnergyRequests energy;
    Box[char[]] others;
}