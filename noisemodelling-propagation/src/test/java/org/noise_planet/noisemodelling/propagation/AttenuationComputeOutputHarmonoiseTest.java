package org.noise_planet.noisemodelling.propagation;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutProfile;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.ProfileBuilder;
import org.noise_planet.noisemodelling.propagation.harmonoise.HarmonoiseAttenuationOutput;
import org.noise_planet.noisemodelling.propagation.harmonoise.HarmonoisePropagationModel;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Tests for Harmonoise prropagation model
 * Ref: Section 3 of Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
 * Acta acustica united with acustica, 97(1), 62-74
 * @author Martin Glesser
 */
public class AttenuationComputeOutputHarmonoiseTest {
    private static final double HUMIDITY = 70;
    private static final double TEMPERATURE = 10;

    public static CutProfile loadCutProfile(String utName) throws IOException {
        //Get test data
        URL url = AttenuationComputeOutputHarmonoiseTest.class.getResource("harmonoise/" + utName + ".json");

        //Out and computation settings
        Assertions.assertNotNull(url);
        try(InputStream inputStream = url.openStream()) {
            ObjectMapper mapper = new ObjectMapper();
            return mapper.readValue(inputStream, CutProfile.class);
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    private static AttenuationOutput computeHarmonoiseAttenuation(CutProfile cutProfile, double radius) {
        //Create profile builder
        ProfileBuilder profileBuilder = new ProfileBuilder()
                .finishFeeding();

        //Set third octave band
        List<Integer> freq_lvl = Arrays.stream(ProfileBuilder.DEFAULT_FREQUENCIES_THIRD_OCTAVE).boxed().collect(Collectors.toList());
        profileBuilder.setFrequencyArray(freq_lvl);

        //Propagation data building
        SceneWithAttenuation sceneWithAttenuation = new SceneWithAttenuation(profileBuilder);
        sceneWithAttenuation.sourceGs.put(-1L, 0.5);

        //Propagation process path data building
        sceneWithAttenuation.defaultCnossosParameters.setHumidity(HUMIDITY);
        sceneWithAttenuation.defaultCnossosParameters.setTemperature(TEMPERATURE);
        sceneWithAttenuation.defaultCnossosParameters.setRadius(radius);

        PropagationModel propagationModel = new HarmonoisePropagationModel();

        return propagationModel.computeAttenuation(sceneWithAttenuation, cutProfile,
                sceneWithAttenuation.defaultCnossosParameters, false).getFirst();
    }

    /**
     * Test case 1 from Harmonoise publication (hard/rigid ground and
     * non refracting atmosphere)
     * Not tested:
     * - diffraction attenuation
     * - Fresnel weighting
     * - modified Fresnel weighting
     */
    @Test
    public void testHarmonoiseCase01() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_01");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                6.036446469248293	,
                6.10478359908884	,
                6.10478359908884	,
                6.10478359908884	,
                6.10478359908884	,
                6.10478359908884	,
                6.10478359908884	,
                6.036446469248293	,
                5.9681093394077465	,
                5.8997722095672	,
                5.831435079726653	,
                5.694760820045559	,
                5.4897494305239185	,
                5.148063781321184	,
                4.533029612756266	,
                3.4396355353075165	,
                1.662870159453302	,
                -1.8906605922551236	,
                -10.979498861047835	,
                -4.760820045558086	,
                2.619589977220958	,
                5.626423690205012	,
                4.3963553530751724	,
                -3.8724373576309787	,
                4.191343963553532	,
                1.457858769931665	,
                4.259681093394079
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.6);
        }
    }

    /**
     * Test case 2 from Harmonoise publication (100 kPa s/m2 ground flow res. and
     * non refracting atmosphere)
     */
    @Test
    public void testHarmonoiseCase02() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_02");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                5.9122401847575	,
                5.77367205542725	,
                5.49653579676674	,
                5.21939953810623	,
                4.80369515011547	,
                4.04157043879907	,
                3.14087759815242	,
                1.82448036951501	,
                0.0230946882217076	,
                -2.40184757505774	,
                -4.89607390300231	,
                -7.04387990762125	,
                -8.15242494226328	,
                -6.76674364896074	,
                -3.92609699769053	,
                -1.08545034642032	,
                1.33949191685912	,
                3.21016166281755	,
                4.59584295612009	,
                4.73441108545034	,
                2.86374133949192	,
                -3.37182448036952	,
                -0.669745958429562	,
                4.38799076212471	,
                0.300230946882216	,
                3.27944572748268	,
                0.577367205542725	,
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.7);
        }
    }

    /**
     * Test case 3 from Harmonoise publication (hard/rigid ground and
     * linear wind profile)
     */
    @Test
    public void testHarmonoiseCase03() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_01");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 340/0.2);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                6.61290322580646	,
                6.61290322580646	,
                6.75115207373272	,
                6.75115207373272	,
                6.75115207373272	,
                6.75115207373272	,
                6.68202764976959	,
                6.61290322580646	,
                6.54377880184332	,
                6.54377880184332	,
                6.40552995391705	,
                6.12903225806452	,
                5.64516129032258	,
                4.74654377880185	,
                3.91705069124424	,
                1.42857142857143	,
                -2.64976958525346	,
                -12.8110599078341	,
                -3.06451612903226	,
                3.64055299539171	,
                6.40552995391705	,
                4.33179723502304	,
                -2.85714285714285	,
                5.43778801843318	,
                0.806451612903228	,
                4.33179723502304	,
                3.64055299539171
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.6);
        }
    }

    /**
     * Test case 5 from Harmonoise publication (2000 kPa.s/m2 ground and
     * non refracting atmosphere)
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 3)
     */
    @Test
    public void testHarmonoiseCase05() throws IOException {
        CutProfile cutProfile = loadCutProfile("case_05");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0);
    }

    /**
     * Test case 7 from Harmonoise publication (200 kPa.s/m2 ground and
     * non refracting atmosphere)
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 3)
     */
    @Test
    public void testHarmonoiseCase07() throws IOException {
        CutProfile cutProfile = loadCutProfile("case_07");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0);
    }

    /**
     * Test case 10 from Harmonoise publication (100 kPa.s/m2 ground and
     * non refracting atmosphere, long-distance propagation)
     */
    @Test
    public void testHarmonoiseCase10() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_09");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                6.199095022624434	,
                6.0633484162895925	,
                5.723981900452486	,
                5.316742081447963	,
                4.638009049773753	,
                3.348416289592759	,
                1.108597285067873	,
                -2.760180995475114	,
                -9.140271493212671	,
                -17.692307692307693	,
                -24.20814479638009	,
                -26.787330316742082	,
                -27.058823529411764	,
                -25.497737556561084	,
                -22.511312217194572	,
                -19.38914027149321	,
                -16.67420814479638	,
                -14.095022624434389	,
                -11.719457013574662	,
                -9.276018099547512	,
                -7.171945701357467	,
                -5.067873303167422	,
                -3.167420814479641	,
                -1.1990950226244372	,
                0.7013574660633459	,
                2.3981900452488674	,
                3.959276018099544
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.3);
        }
    }

    /**
     * Test case 12 from Harmonoise publication (rigid ground and
     * non refracting atmosphere, long-distance propagation)
     */
    @Test
    public void testHarmonoiseCase12() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_11");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                6.004514672686234	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.004514672686234	,
                6.0722347629796865	,
                6.0722347629796865	,
                6.0722347629796865	,
                5.869074492099326	,
                5.936794582392778	,
                5.869074492099326	,
                5.665914221218962	,
                5.665914221218962	,
                5.59819413092551	,
                5.191873589164789	,
                4.717832957110613	,
                4.040632054176076	,
                2.889390519187362	,
                1.0609480812641117	,
                -1.9187358916478523
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.6);
        }
    }
}
