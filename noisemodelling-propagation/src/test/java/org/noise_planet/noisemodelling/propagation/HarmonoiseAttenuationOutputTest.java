package org.noise_planet.noisemodelling.propagation;

import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;
import org.noise_planet.noisemodelling.propagation.harmonoise.HarmonoiseAttenuation;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HarmonoiseAttenuationOutputTest {
    /**
     * Test the computation of the boundary loss factor in grazing propagation
     * Ref: K. Attenborough, K. M. Li, and K. Horoshenkov, Predicting Outdoor Sound. Taylor & Francis, 2006.
     * doi: 10.1201/9781482295023. (figure 2.6)
     *
     */
    @Test
    public void boundaryLossFactorTest(){
        double[] magnitude = new double[]{0.001, 0.002, 0.004, 0.006, 0.01, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6,
                1, 2, 4, 6, 10, 20};
        double[] phaseAngle = new double[]{0, Math.PI/4, 80 * Math.PI / 180, Math.PI/2};
        double[][] boundaryLossFactor = new double[phaseAngle.length][magnitude.length];
        for (int i = 0; i < phaseAngle.length; i++){
            int finalI = i;
            boundaryLossFactor[i] = Arrays.stream(magnitude)
                    .map(a -> HarmonoiseAttenuation.boundaryLossFactor(
                            new Complex(0, phase(phaseAngle[finalI])).exp().multiply(a)
                    ).abs())
                    .toArray();
        }
        double[] expected_0 = new double[]{
                0.014641288433381305	,
                0.2196193265007338	,
                -0.08784773060029494	,
                -0.2928257686676403	,
                -0.7027818448023382	,
                -2.035139092240115	,
                -4.187408491947288	,
                -6.032210834553439	,
                -10.131771595900439	,
                -18.84333821376281	,
                -30.01464128843338	,
                -36.67642752562225	,
                -45.900439238653
        }; // plotdigitized from Attenborough et al., distance from 0.01 to 10
        for (int i = 0; i < expected_0.length; i++) {
            assertEquals(expected_0[i], 20 * Math.log10( boundaryLossFactor[0][i+4]), 1.5);
        }
        double[] expected_90 = new double[]{
                0.19061583577713037	,
                0.395894428152495	,
                0.4985337243401773	,
                0.9090909090909136	,
                1.3196480938416428	,
                2.34604105571848	,
                4.604105571847512	,
                6.6568914956011795	,
                10.659824046920825	,
                17.33137829912023
        }; // plotdigitized from Attenborough et al., distance from 0.01 to 2
        for (int i = 0; i < expected_90.length; i++) {
            assertEquals(expected_90[i], 20 * Math.log10( boundaryLossFactor[3][i+4]), 1);
        }
    }

    /**
     * Return the phase of the numerical distance for a given phase angle of impedance in the case
     * of grazing propagation
     * Ref: K. Attenborough, K. M. Li, and K. Horoshenkov, Predicting Outdoor Sound. Taylor & Francis, 2006.
     * doi: 10.1201/9781482295023. (equation 2.68b)
     *
     * @param angle phase angle of surface impedance [rad]
     * @return phase of the numerical distance
     */
    private double phase(double angle){
        return Math.atan( Math.cos(Math.PI/4 + angle) / Math.cos(Math.PI/4 - angle) );
    }
}
