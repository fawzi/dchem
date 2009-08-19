module dchem.sys.ReadInStruct;

struct Particle{
    PIndex pIndex;
    PIndex molIndex;
    char[] name;
    char[] molName;
    real[3] pos;
}

convertToSys
