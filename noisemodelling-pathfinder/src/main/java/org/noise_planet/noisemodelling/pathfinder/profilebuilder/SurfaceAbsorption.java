/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */
package org.noise_planet.noisemodelling.pathfinder.profilebuilder;

import org.noise_planet.noisemodelling.pathfinder.utils.ComplexNumber;

/**
 * Collection of methods related to wall and ground absorption coefficients
 * @author Nicolas Fortin
 * @author Martin Glesser
 */
public class SurfaceAbsorption {

    /**
     * Get WallAlpha
     */
    public static double getWallAlpha(double wallAlpha, double freq_lvl)
    {
        double value;
        if(wallAlpha >= 0 && wallAlpha <= 1) {
            // todo let the user choose if he wants to convert G to Sigma
            //value = GetWallImpedance(20000 * Math.pow (10., -2 * Math.pow (wallAlpha, 3./5.)),freq_lvl);
            value= wallAlpha;
        } else {
            value = computeSurfaceImpedance(Math.min(20000, Math.max(20, wallAlpha)),freq_lvl);
        }
        return value;
    }

    /**
     * Compute the surface normalized acoustic impedance using Delany-Bazley ground impedance models
     * Ref: Delany M. E. and Bazley E. N., Acoustical properties of fibrous absorbent materials,
     * Applied Acoustics 3, 1970, pp. 105-116
     *
     * @param sigma surface flow resistivity [kPa.s/m2]
     * @param freq_l frequency [Hz]
     * @return normalized ground impedance[]
     */
    public static double computeSurfaceImpedance(double sigma, double freq_l){
        return computeSurfaceImpedance(sigma, freq_l, "delany-bazley");
    }

    /**
     * Compute the surface normalized acoustic impedance using different ground impedance models
     *
     * @param sigma surface flow resistivity [kPa.s/m2]
     * @param freq_l frequency [Hz]
     * @param method ground impedance model
     * @return normalized ground impedance[]
     */
    public static double computeSurfaceImpedance(double sigma, double freq_l, String method)
    {
        ComplexNumber Z;
        switch (method) {
            default:
                // Ref: Delany M. E. and Bazley E. N., Acoustical properties of fibrous absorbent
                // materials, Applied Acoustics 3, 1970, pp. 105-116
                double s = Math.log(freq_l / sigma);
                double x = 1. + 9.08 * Math.exp(-.75 * s);
                double y = - 11.9 * Math.exp(-0.73 * s);
                Z = new ComplexNumber(x, y);
                break;
                // Another fit giving similar results is sometimes found in the literature with:
                // x = 1 + 0.0571 * (rho_0 * s)^-0.754
                // y = 0.0870 * (rho_0 * s)^-0.732
        }
        return computeGroundAbsorptionCoefficient(Z);
    }

    /**
     * Compute the Ground absorption coefficient
     *
     * @param impedance normalized ground impedance[]
     * @return ground absorption coefficient []
     */
    static double computeGroundAbsorptionCoefficient(ComplexNumber impedance)         // TODO convert impedance to alpha
    {
        double alpha ;
        ComplexNumber z = ComplexNumber.divide(new ComplexNumber(1.0,0), impedance) ;
        double x = z.getRe();
        double y = z.getIm();
        double a1 = (x * x - y * y) / y ;
        double a2 = y / (x * x + y * y + x) ;
        double a3 = ((x + 1) *(x + 1) + y * y) / (x * x + y * y) ;
        alpha = 8 * x * (1 + a1 * Math.atan(a2) - x * Math.log(a3)) ;
        return alpha ;
    }
}
