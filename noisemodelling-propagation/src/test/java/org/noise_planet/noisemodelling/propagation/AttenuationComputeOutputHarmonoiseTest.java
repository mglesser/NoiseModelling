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
import java.nio.file.Paths;
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
    double[] DEFAULT_FREQUENCIES_THIRD_OCTAVE = new double[] {50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000};
    String rootDir = Paths.get(System.getProperty("user.dir")).getParent().toString();

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

    private static AttenuationOutput computeHarmonoiseAttenuation(CutProfile cutProfile, double aMeteo, double bMeteo) {
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
        sceneWithAttenuation.defaultCnossosParameters.setaMeteo(aMeteo);
        sceneWithAttenuation.defaultCnossosParameters.setbMeteo(bMeteo);
        sceneWithAttenuation.defaultCnossosParameters.setZ0(0.1);

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
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 0);
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
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case01.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.7);
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
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 0);
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
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case02.svg");
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
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0.2, 0);
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
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case03.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 4);
        }
    }

    /**
     * Test case 4 from Harmonoise publication (100 kPa s/m2 ground flow res. ground and
     * linear celerity profile)
     */
    @Test
    public void testHarmonoiseCase04() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_02");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0.2, 0);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                6.19815668202765	,
                6.19815668202765	,
                6.05990783410138	,
                5.36866359447005	,
                5.78341013824885	,
                4.60829493087558	,
                3.57142857142857	,
                2.1889400921659	,
                0.046082949308758	,
                -2.30414746543779	,
                -4.72350230414746	,
                -7.14285714285714	,
                -7.69585253456221	,
                -5.4147465437788	,
                -1.82027649769585	,
                1.29032258064516	,
                3.50230414746544	,
                4.67741935483871	,
                4.12442396313364	,
                1.5668202764977	,
                -4.10138248847926	,
                0.875576036866363	,
                4.19354838709678	,
                -0.161290322580644	,
                3.70967741935484	,
                1.3594470046083	,
                2.1889400921659
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case04.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 3);
        }
    }

    /**
     * Test case 5 from Harmonoise publication (Sound barrrier, 20000 kPa.s/m2 ground and
     * non refracting atmosphere)
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 3)
     */
    @Test
    public void testHarmonoiseCase05() throws IOException {
        CutProfile cutProfile = loadCutProfile("case_05");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 0);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                1.9529411764705884	,
                1.0352941176470587	,
                0.2588235294117638	,
                -0.8705882352941181	,
                -2.5647058823529427	,
                -4.8235294117647065	,
                -8.070588235294117	,
                -12.658823529411766	,
                -13.505882352941176	,
                -9.83529411764706	,
                -8.988235294117647	,
                -12.023529411764706	,
                -22.823529411764707	,
                -24.870588235294118	,
                -17.24705882352941	,
                -13.647058823529413	,
                -14.63529411764706	,
                -18.870588235294118	,
                -20.705882352941178	,
                -17.88235294117647	,
                -22.68235294117647	,
                -21.83529411764706	,
                -22.611764705882354	,
                -23.95294117647059	,
                -25.011764705882353	,
                -26.28235294117647	,
                -27.411764705882355
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case05.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 3);
        }
    }

    /**
     * Test case 7 from Harmonoise publication (200 kPa.s/m2 ground and
     * non refracting atmosphere)
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 3)
     */
//    @Test
//    public void testHarmonoiseCase07() throws IOException {
//        CutProfile cutProfile = loadCutProfile("case_07");
//        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 0);
//        assert(output.aGlobal.length !=0);
//    }

    /**
     * Test case 9 from Harmonoise publication (100 kPa.s/m2 ground and
     * linear celerity profile, long-distance propagation)
     */
    @Test
    public void testHarmonoiseCase09() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_09");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 1);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                9.458239277652371	,
                9.59367945823928	,
                9.729119638826187	,
                9.458239277652371	,
                8.713318284424382	,
                6.410835214446955	,
                3.566591422121899	,
                -0.9029345372460469	,
                -6.523702031602708	,
                -10.383747178329571	,
                -10.383747178329571	,
                -7.268623024830699	,
                -4.085778781038373	,
                0.38374717832957117	,
                4.717832957110613	,
                6.275395033860047	,
                2.2799097065462774	,
                -1.3092550790067712	,
                7.291196388261852	,
                7.832957110609485	,
                3.1602708803611748	,
                -1.6478555304740397	,
                5.462753950338602	,
                7.155756207674944	,
                2.0090293453724613	,
                7.2234762979684	,
                4.853273137697521
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case09.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 1);
        }
    }

    /**
     * Test case 10 from Harmonoise publication (100 kPa.s/m2 ground and
     * non refracting atmosphere, long-distance propagation)
     */
    @Test
    public void testHarmonoiseCase10() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_09");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 0);
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
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case10.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.3);
        }
    }

    /**
     * Test case 11 from Harmonoise publication (rigid ground and
     * linear celerity profile, long-distance propagation)
     */
    @Test
    public void testHarmonoiseCase11() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_11");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 1);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                9.932126696832576	,
                10.746606334841626	,
                11.357466063348415	,
                12.239819004524882	,
                12.850678733031671	,
                13.529411764705877	,
                14.004524886877824	,
                14.004524886877824	,
                14.072398190045249	,
                13.665158371040718	,
                13.257918552036195	,
                12.036199095022624	,
                9.117647058823525	,
                3.348416289592759	,
                -2.4886877828054317	,
                1.3122171945701346	,
                8.167420814479637	,
                7.76018099547511	,
                -8.597285067873305	,
                2.8733031674208114	,
                9.25339366515837	,
                10.610859728506785	,
                7.69230769230769	,
                4.434389140271492	,
                10.135746606334838	,
                5.384615384615383	,
                8.50678733031674
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case11.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 1);
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
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 0);
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
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case12.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.2);
        }
    }

    /**
     * Test case 14 from Harmonoise publication (Sound barrier, 100 kPa.s/m2 ground
     * non refracting atmosphere and long-distance propagation)
     */
    @Test
    public void testHarmonoiseCase14() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_13");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 0);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                2.589928057553955	,
                1.7985611510791344	,
                0.7194244604316538	,
                -0.6474820143884905	,
                -2.6618705035971217	,
                -5.539568345323742	,
                -9.20863309352518	,
                -14.244604316546763	,
                -21.07913669064748	,
                -25.39568345323741	,
                -27.841726618705035	,
                -28.633093525179856	,
                -26.690647482014388	,
                -22.949640287769785	,
                -18.992805755395683	,
                -16.618705035971225	,
                -16.546762589928058	,
                -19.35251798561151	,
                -16.618705035971225	,
                -15.899280575539569	,
                -20.64748201438849	,
                -23.165467625899282	,
                -23	,
                -22.66187050359712	,
                -22.949640287769785	,
                -28.561151079136692	,
                -25.323741007194243
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case14.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 0.2);
        }
    }

    /**
     * Test case 16 from Harmonoise publication (Sound barrier, rigid ground
     * non refracting atmosphere and long-distance propagation)
     */
    @Test
    public void testHarmonoiseCase16() throws IOException {
        // Load data and compute attenuation
        CutProfile cutProfile = loadCutProfile("case_15");
        AttenuationOutput output = computeHarmonoiseAttenuation(cutProfile, 0, 0);
        assert (output instanceof HarmonoiseAttenuationOutput);
        //Assertion
        double[] referenceExcessAttenuation = {
                3.221153846153843	,
                2.7163461538461533	,
                1.9230769230769198	,
                1.41826923076923	,
                0.7692307692307665	,
                -0.02403846153846345	,
                -0.7451923076923102	,
                -1.3221153846153868	,
                -2.76442307692308	,
                -3.8461538461538467	,
                -5.721153846153847	,
                -8.245192307692308	,
                -12.139423076923077	,
                -18.629807692307693	,
                -21.58653846153846	,
                -15.024038461538462	,
                -13.653846153846153	,
                -16.53846153846154	,
                -29.73557692307692	,
                -31.3	,
                -21.08173076923077	,
                -18.91826923076923	,
                -19.927884615384617	,
                -24.326923076923077	,
                -25.552884615384613	,
                -22.8125	,
                -27.139423076923077
        }; // plotdigitized from publication, 1/3 oct band 25Hz - 10kHz
        assert(output.aGlobal.length !=0);
        Graph2D plot = new Graph2D(DEFAULT_FREQUENCIES_THIRD_OCTAVE,
                Arrays.copyOfRange(referenceExcessAttenuation, 3, referenceExcessAttenuation.length),
                output.aGlobal);
        plot.legend("Salomons 2011", "NoiseModelling");
        plot.draw(rootDir + "/Docs/images/Validation/Harmonoise/case16.svg");
        for (int i = 0; i < output.aGlobal.length; i++) {
            assertEquals(referenceExcessAttenuation[i+3], output.aGlobal[i], 1);
        }
    }
}
