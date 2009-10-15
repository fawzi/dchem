/// numbers that guarantee reproducible results at least within the same run
///
/// basic ops +-*/ are available, and it is possible to convert reals and doubles
/// to them. The actual precision is not necessarily high.
///
/// on some platform it could be sostitued by float or double or real if it 
/// has the correct characteristic, by default it uses a fixed point implementation.
///
/// conversion from a different vector type v is simply vec3fR.from(v) (back conversion)
/// is similar: vec3R.from(v)...
///
/// author: fawzi
module dchem.util.Reproducible;
import xf.omg.core.Fixed;
import xf.omg.core.LinearAlgebra;

// useful for conversions to/from fixedReal,fixedReduced (constant, non constant)
public import xf.omg.Algebra: cscalar,scalar;
// usage: cscalar!(type,value) and scalar!(type)(value)

/// type useful for real positions
alias fixedT!(24,8) fixedReal;
/// type useful for reduced positions
alias fixedT!(16,16) fixedReduced;

alias Vector!(fixedReal, 2) vec2fR;
alias Vector!(fixedReal, 3) vec3fR;
alias Vector!(fixedReal, 4) vec4fR;

alias Vector!(fixedReduced, 2)  vec2fr;
alias Vector!(fixedReduced, 3)  vec3fr;
alias Vector!(fixedReduced, 4)  vec4fr;

