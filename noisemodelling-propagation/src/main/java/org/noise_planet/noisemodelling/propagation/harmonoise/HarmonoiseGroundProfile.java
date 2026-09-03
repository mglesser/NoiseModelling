/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.propagation.harmonoise;

import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutPoint;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutProfile;
import org.noise_planet.noisemodelling.pathfinder.utils.ComplexNumber;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.Math.*;

/**
 * 2D ground profile between a source and a receiver.
 * @author Martin Glesser
 */

public class HarmonoiseGroundProfile {
    LineString profile; // 2D coordinates of ground profile vertices
    Coordinate source; // 3D coordinates of the source
    Coordinate receiver; // 3D coordinates of the receiver
    List<Coordinate> vertices = new ArrayList<>();

    /**
     * Initialize HarmonoiseGroundProfile object from CutProfile object.
     *
     * @param cutProfile 3D profile from source to receiver
     */
    public HarmonoiseGroundProfile(CutProfile cutProfile){
        source = cutProfile.getSource().getCoordinate();
        receiver = cutProfile.getReceiver().getCoordinate();
        // Get the whole 2D profile including ground points
        List<Integer> hullIndices = cutProfile.getConvexHullIndices(cutProfile.computePts2D());
        Coordinate[] coordinates = cutProfile.computePts2DGround(hullIndices).toArray(new Coordinate[0]);
        GeometryFactory geometryFactory = new GeometryFactory();
        profile = geometryFactory.createLineString(coordinates);
    }

    /**
     * Generate a curved profile from a coordinate list, two endpoints (source and receiver) and a curvature radius.
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 2.5)
     * Note: This implementation yield similar results to the one from CurvedProfileGenerator.applyTransformation.
     * However, it works only on the whole ground profile (from zGroundSource to zGroundReceiver).
     */
    public void computeCurvedProfile(double radius){
        // Segment Profile (second paragraph of section 2.5)
        double dsr = source.distance(receiver);
        double maxSegmentLength = min( dsr/3 , max(50, dsr/20));
        profile = (LineString) Densifier.densify(profile, maxSegmentLength);

        // Ground curvature
        double hSource = source.z;
        double hReceiver = receiver.z;
        double hm = (hSource + hReceiver) / 2;
        double c0 = 2* (hm + radius); // Eq. 77
        ComplexNumber c = new ComplexNumber(0, c0); // Eq. 76
        double xc = 0.5 * (profile.getStartPoint().getX() + profile.getEndPoint().getX());
        double yc = 0.5 * (profile.getStartPoint().getY() + profile.getEndPoint().getY()) + hm;
        ComplexNumber w0 = new ComplexNumber(xc, yc); // Eq. 75
        double deltaY = 0;
        vertices = Arrays.asList(new Coordinate[profile.getNumPoints()]);
        for (int i = 0; i < profile.getNumPoints(); i++) {
            ComplexNumber w = new ComplexNumber(profile.getCoordinateN(i).getX(), profile.getCoordinateN(i).getY());
            ComplexNumber wPrim = ComplexNumber.divide(
                    ComplexNumber.multiply(c, ComplexNumber.subtract(w, w0)),
                    ComplexNumber.add(c, ComplexNumber.subtract(w, w0))
            ); // Eq. 74

            // Create new coordinate with transformed z (incl. profile translation)
            if (i == 0) {
                deltaY = profile.getCoordinateN(i).getY() - wPrim.getIm();
                vertices.set(i,
                        new Coordinate(wPrim.getRe() + xc, profile.getCoordinateN(i).getY() , profile.getCoordinateN(i).getZ()));
            } else {
                vertices.set(i,
                        new Coordinate(wPrim.getRe() + xc, wPrim.getIm() + deltaY, profile.getCoordinateN(i).getZ()));
            }
        }
    }
}
