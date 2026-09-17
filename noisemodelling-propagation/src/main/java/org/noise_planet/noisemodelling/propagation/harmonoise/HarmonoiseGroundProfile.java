/**
 * NoiseModelling is a library capable of producing noise maps. It can be freely used either for research and education, as well as by experts in a professional use.
 * <p>
 * NoiseModelling is distributed under GPL 3 license. You can read a copy of this License in the file LICENCE provided with this software.
 * <p>
 * Official webpage : http://noise-planet.org/noisemodelling.html
 * Contact: contact@noise-planet.org
 */

package org.noise_planet.noisemodelling.propagation.harmonoise;

import org.apache.commons.math3.complex.Complex;
import org.locationtech.jts.algorithm.Angle;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.CutProfile;

import java.util.Arrays;
import java.util.List;

import static java.lang.Math.*;

/**
 * 2D ground profile between a source and a receiver.
 * @author Martin Glesser
 */

public class HarmonoiseGroundProfile {
    private final Coordinate[] vertices;

    /**
     * Initialize HarmonoiseGroundProfile object from CutProfile object.
     *
     * @param cutProfile 3D profile from source to receiver
     */
    public HarmonoiseGroundProfile(CutProfile cutProfile){

        this(cutProfile, 0);
    }

    /**
     * Initialize HarmonoiseGroundProfile object from CutProfile object.
     *
     * @param cutProfile 3D profile from source to receiver
     * @param radius curvature radius of the equivalent ground profile
     */
    public HarmonoiseGroundProfile(CutProfile cutProfile, double radius){
        // Get the whole 2D profile including ground points
        List<Integer> hullIndices = cutProfile.getConvexHullIndices(cutProfile.computePts2D());
        Coordinate[] groundVertices = cutProfile.computePts2DGround(hullIndices).toArray(new Coordinate[0]);
        if (radius != 0){
            vertices = computeCurvedProfile(groundVertices, radius);
        } else {
            vertices = groundVertices;
        }
        vertices[0].y = cutProfile.getSource().getCoordinate().getZ();
        vertices[vertices.length-1].y = cutProfile.getReceiver().getCoordinate().getZ();
    }

    /**
     * Generate a curved profile from a coordinate list, two endpoints (source and receiver) and a curvature radius.
     * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
     * Acta acustica united with acustica, 97(1), 62-74 (section 2.5)
     * Note: This implementation yield similar results to the one from CurvedProfileGenerator.applyTransformation.
     * However, it works only on the whole ground profile (from zGroundSource to zGroundReceiver).
     */
    private Coordinate[] computeCurvedProfile(Coordinate[] groundProfile, double radius){
        GeometryFactory geometryFactory = new GeometryFactory();
        LineString profile = geometryFactory.createLineString(groundProfile);
        Coordinate source = groundProfile[0];
        Coordinate receiver = groundProfile[groundProfile.length-1];
        Coordinate[] vertices;
        // Segment Profile (second paragraph of section 2.5)
        double dsr = source.distance(receiver);
        double maxSegmentLength = min( dsr/3 , max(50, dsr/20));
        profile = (LineString) Densifier.densify(profile, maxSegmentLength);

        // Ground curvature
        double hSource = source.z;
        double hReceiver = receiver.z;
        double hm = (hSource + hReceiver) / 2;
        double c0 = 2* (hm + radius); // Eq. 77
        Complex c = new Complex(0, c0); // Eq. 76
        double xc = 0.5 * (profile.getStartPoint().getX() + profile.getEndPoint().getX());
        double yc = 0.5 * (profile.getStartPoint().getY() + profile.getEndPoint().getY()) + hm;
        Complex w0 = new Complex(xc, yc); // Eq. 75
        double deltaY = 0;
        vertices = new Coordinate[profile.getNumPoints()];
        for (int i = 0; i < profile.getNumPoints(); i++) {
            Complex w = new Complex(profile.getCoordinateN(i).getX(), profile.getCoordinateN(i).getY());
            Complex wPrim = c.multiply(w.subtract(w0)).divide(c.add(w.subtract(w0))); // Eq. 74

            // Create new coordinate with transformed z (incl. profile translation)
            if (i == 0) {
                deltaY = profile.getCoordinateN(i).getY() - wPrim.getImaginary();
                vertices[i] =
                        new Coordinate(wPrim.getReal() + xc, profile.getCoordinateN(i).getY() , profile.getCoordinateN(i).getZ());
            } else {
                vertices[i] =
                        new Coordinate(wPrim.getReal() + xc, wPrim.getImaginary() + deltaY, profile.getCoordinateN(i).getZ());
            }
        }
        return vertices;
    }

    public Coordinate[] getVertices() {
        return vertices;
    }

    public int getNVertices() {
        return vertices.length;
    }

    /**
     * @return ground profile vertices with first and last ground vertices replaced
     * respectively by the source and the receiver vertices
     */
    public Coordinate[] getVertices(int iStart, int iEnd) {
        return Arrays.copyOfRange(vertices, iStart, iEnd+1);
    }

    /**
     * Getter for vertex coordinate.
     *
     * @param i index of the vertex
     * @return coordinate of the vertex
     */
    public Coordinate getVertex(int i){
        return vertices[i];
    }

    /**
     * Return the coordinates of the image of a vertex with respect to a ground segment plane
     *
     * @param iVertex index of the vertex
     * @param iSeg index of the first point of the ground segment
     * @return image of the vertex
     */
    public Coordinate getImageVertex(int iVertex, int iSeg) {
        double x0 = vertices[iVertex].getX();
        double y0 = vertices[iVertex].getY();
        double x1 = vertices[iSeg].getX();
        double x2 = vertices[iSeg+1].getX();
        double y1 = vertices[iSeg].getY();
        double y2 = vertices[iSeg+1].getY();
        // Segment line equation: ax + by + c = 0
        double a = 1 / (x2 - x1);
        double b = 1 / (y1 - y2);
        double c = y1 / (y2 - y1) - x1 / (x2 - x1);
        // Image vertex coordinates
        double xi = x0 - 2*a * (a*x0 + b*y0 + c) / (a*a + b*b);
        double yi = y0 - 2*b * (a*x0 + b*y0 + c) / (a*a + b*b);
        return new Coordinate(xi, yi);
    }

    /**
     * Return the height of the (secondary) source relative to a ground segment plane
     * @param iSeg index of the ground segment
     * @param iSrc index of the (secondary) source
     * @return local height of the source
     */
    public double getLocalSourceHeight(int iSeg, int iSrc){
        Coordinate source = vertices[iSrc];
        Coordinate segmentEnd = vertices[iSeg+1];
        return segmentEnd.distance(source)
                * Math.sin(Angle.angleBetweenOriented(source, segmentEnd, vertices[iSeg]));
    }

    /**
     * Return the height of the (secondary) receiver relative to a ground segment plane
     * @param iSeg index of the ground segment
     * @param iRcv index of the (secondary) receiver
     * @return local height of the receiver
     */
    public double getLocalReceiverHeight(int iSeg, int iRcv) {
        Coordinate receiver = vertices[iRcv];
        Coordinate segmentStart = vertices[iSeg];
        return segmentStart.distance(receiver)
                * Math.sin(Angle.angleBetweenOriented(vertices[iSeg + 1], segmentStart, receiver));
    }
}
