package org.noise_planet.noisemodelling.propagation.harmonoise;


import org.noise_planet.noisemodelling.propagation.PropagationModel;
import org.noise_planet.noisemodelling.propagation.PropagationModelFactory;
import org.noise_planet.noisemodelling.propagation.template.TemplatePropagationModel;

/**
 * Declares the concrete factory method that returns TemplatePropagationModel objects
 * @author Martin Glesser
 */
public class HarmonoisePropagationModelFactory implements PropagationModelFactory {
    /**
     * Factory method that returns TemplatePropagationModel objects
     * @return PropagationModel object
     */
    public PropagationModel create(){
        return new TemplatePropagationModel();
    }
}
