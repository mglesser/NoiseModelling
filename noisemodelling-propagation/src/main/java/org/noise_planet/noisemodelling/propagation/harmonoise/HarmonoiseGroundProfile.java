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
import org.locationtech.jts.algorithm.Intersection;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.math.Vector2D;
import org.noise_planet.noisemodelling.pathfinder.profilebuilder.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.Math.*;

/**
 * 2D ground profile between a source and a receiver.
 * Ref: Salomons, E., Van Maercke, D., Defrance, J.,&amp;De Roo, F. (2011). The Harmonoise sound propagation model.
 * Acta acustica united with acustica, 97(1), 62-74
 * @author Martin Glesser
 */

public class HarmonoiseGroundProfile {
    private Coordinate[] vertices;
    private double[] flowResistivity;
    private Complex[][] reflectionCoefficient;
    private double sourceHeight;
    private double receiverHeight;

    /**
     * Initialize HarmonoiseGroundProfile object from CutProfile object.
     *
     * @param cutProfile 3D profile from source to receiver
     * @param nFreq size of the frequency axis
     */
    public HarmonoiseGroundProfile(CutProfile cutProfile, int nFreq){this(cutProfile, nFreq,
            0, 0, 0, 0);
    }

    /**
     * Initialize HarmonoiseGroundProfile object from CutProfile object considering
     * log/lin celerity profile: c(z) = c0 + az + bln(1+z/z0).
     *
     * @param cutProfile 3D profile from source to receiver
     * @param nFreq size of the frequency axis
     * @param aMeteo linear component of the lin/log celerity profile
     * @param bMeteo log component of the lin/log celerity profile
     */
    public HarmonoiseGroundProfile(CutProfile cutProfile, int nFreq, double aMeteo, double bMeteo,
                                   double z0, double c0){
        // Get the whole 2D profile including ground points
        extractGroundVertices(cutProfile);
        double radius = computeRadius(aMeteo, bMeteo, z0, c0);
        if(radius != 0) {
            double dsr = cutProfile.getCutPoints().getFirst().getCoordinate()
                    .distance(cutProfile.getCutPoints().getLast().getCoordinate());
            refine(dsr);
            curve(radius);
        }
        // The real source and real receiver ground points are replaced by respectively, the source and the receiver
        vertices[0].setY(vertices[0].getY() + sourceHeight);
        vertices[getNVertices() - 1].setY(vertices[getNVertices() - 1].getY() + receiverHeight);
        reflectionCoefficient = new Complex[getNVertices()-1][nFreq];
    }

    /**
     * Extract 2D vertices and segments flow resistivity from the cutProfile
     *
     * @param cutProfile 3D profile from source to receiver
     */
    private void extractGroundVertices(CutProfile cutProfile){
        // Extract top-elevation point, including ground points,
        List<Integer> indices = new ArrayList<>(0);
        vertices = cutProfile.computePts2DGround(0, indices).toArray(new Coordinate[0]);
        int nVertices = vertices.length;
        // source and receiver height from ground
        sourceHeight = cutProfile.getCutPoints().getFirst().getCoordinate().getZ() - vertices[0].getY();
        receiverHeight = cutProfile.getCutPoints().getLast().getCoordinate().getZ() - vertices[nVertices-1].getY();
        // Extract flow resistivity
        flowResistivity = extractFlowResistivity(cutProfile, nVertices, indices);
    }

    /**
     * Extract ground segment's flow resistivity from cut profile
     *
     * @param cutProfile 3D profile from source to receiver
     * @param nVertices ground profile's number of vertices
     * @param indices correspondence between cut points and vertices indices
     * @return ground segment's flow resistivity
     */
    private double[] extractFlowResistivity(CutProfile cutProfile, int nVertices, List<Integer> indices){
        double[] resistivity = new double[nVertices];
        for (int i = 0; i < indices.size(); i++) {
            if (cutProfile.cutPoints.get(i) instanceof CutPointWall) {
                // Set rigid wall
                // TODO : get absorption from cutProfile
                resistivity[indices.get(i)-2] = 1e6;
                resistivity[indices.get(i)-1] = 1e6;
                resistivity[indices.get(i)] = cutProfile.cutPoints.get(i).groundCoefficient;
            } else {
                resistivity[indices.get(i)] = cutProfile.cutPoints.get(i).groundCoefficient;
            }
        }
        return resistivity;
    }

    /**
     * Refine the segmentation of the profile according to section 2.5 of Salomons et al.
     *
     * @param dsr distance source-receiver
     */
    private void refine(double dsr){
        GeometryFactory geometryFactory = new GeometryFactory();
        double maxSegmentLength = min( dsr/3 , max(50, dsr/20));
        Coordinate[] refinedVertices = new Coordinate[] {vertices[0]};
        double[] refinedResistivity = new double[] {flowResistivity[0]};
        for (int k = 0; k < getNVertices() - 1; k++) {
            // Refine current segment
            LineString segment = geometryFactory.createLineString(
                    new Coordinate[]{vertices[k], vertices[k + 1]});
            segment = (LineString) Densifier.densify(segment, maxSegmentLength);
            Coordinate[] verticesConcat = new Coordinate[refinedVertices.length + segment.getCoordinates().length - 1];
            System.arraycopy(refinedVertices, 0, verticesConcat, 0, refinedVertices.length);
            System.arraycopy(Arrays.copyOfRange(segment.getCoordinates(), 1, segment.getCoordinates().length),
                    0, verticesConcat, refinedVertices.length, segment.getCoordinates().length - 1);
            refinedVertices = verticesConcat;
            // Update flow resistivity vector accordingly
            double[] resistivity = new double[segment.getCoordinates().length-1];
            Arrays.fill(resistivity, flowResistivity[k]);
            double[] resistivityConcat = new double[refinedResistivity.length + resistivity.length];
            System.arraycopy(refinedResistivity, 0, resistivityConcat, 0, refinedResistivity.length);
            System.arraycopy(resistivity,0, resistivityConcat, refinedResistivity.length, resistivity.length);
            refinedResistivity = resistivityConcat;
        }
        vertices = refinedVertices;
        flowResistivity = refinedResistivity;
    }

    /**
     * Compute the equivalent curvature radius Rc of the linear profile c(z) = c0 + c0/Rc.z)
     * that agree as much as possible with the lin/log profile c(z) = c0 + az + bln(1+z/z0)
     * according to Harmonoise methodology.
     * Ref: Salomons et al section 3
     *
     * @return the equivalent curvature radius
     */
    public double computeRadius(double aMeteo, double bMeteo, double z0, double c0){
        if (aMeteo ==0 && bMeteo ==0) { return 0;}
        double hm = (sourceHeight + receiverHeight) / 2;
        // Linear component of the celerity profile
        double ra = aMeteo / c0;
        // Logarithmic component of the celerity profile
        double rb;
        if(bMeteo != 0){
            if (bMeteo > 0) {
                double vDist = sourceHeight - receiverHeight;
                double hDist = vertices[vertices.length-1].getX() - vertices[0].getX();
                double tan = vDist / hDist ;
                double distance = Math.sqrt(vDist*vDist + hDist*hDist);
                double k = Math.sqrt(bMeteo / (2 * Math.PI * c0));
                double gamma = (1 + 4 * k * k)/(1 - 4 * k * k);
                double a = 1 + tan * tan - gamma * gamma;
                double b = hm * (1 + tan * tan);
                double c = b * hm + distance * distance / 4;
                double d = b * b - a * c;
                if (d > 0){
                    rb = (Math.sqrt(d) - b) / c ;
                    if (rb < 0){ rb = 0; }
                }
                else {
                    rb = 0;
                }
            }
            else {
                rb = 0;
            }
        }
        else {
            rb = bMeteo / c0 / hm ;
        }
        return Math.max(1 / (ra + rb), 1.0e-6 );
    }

    /**
     * Generate a curved profile from a coordinate list, two endpoints (source and receiver) and a curvature radius.
     * Ref: Salomons et al. (section 2.5)
     * Note: This implementation yield similar results to the one from CurvedProfileGenerator.applyTransformation.
     * However, it works only on the whole ground profile (from zGroundSource to zGroundReceiver).
     */
    private void curve(double radius){
        int nVertices = vertices.length;
        Coordinate sourceGroundPoint = vertices[0];
        Coordinate receiverGroundPoint = vertices[nVertices - 1];
        Coordinate[] curvedVertices;
        // Ground curvature
        double hm = (sourceHeight + receiverHeight) / 2;
        double c0 = 2 * (hm + radius); // Eq. 77
        Complex c = new Complex(0, c0); // Eq. 76
        double xc = 0.5 * (sourceGroundPoint.getX() + receiverGroundPoint.getX());
        double yc = 0.5 * (sourceGroundPoint.getY() + receiverGroundPoint.getY()) + hm;
        Complex w0 = new Complex(xc, yc); // Eq. 75
        double deltaY = 0;
        curvedVertices = new Coordinate[nVertices];
        for (int i = 0; i < nVertices; i++) {
            Complex w = new Complex(vertices[i].getX(), vertices[i].getY());
            Complex wPrim = c.multiply(w.subtract(w0)).divide(c.add(w.subtract(w0))); // Eq. 74
            // Create new coordinate with transformed z (incl. profile translation)
            if (i == 0) {
                deltaY = sourceGroundPoint.getY() - wPrim.getImaginary();
                curvedVertices[i] =
                        new Coordinate(wPrim.getReal() + xc, vertices[i].getY() , vertices[i].getZ());
            } else {
                curvedVertices[i] =
                        new Coordinate(wPrim.getReal() + xc, wPrim.getImaginary() + deltaY, vertices[i].getZ());
            }
        }
        vertices = curvedVertices;
    }

    /**
     * @return number of vertices
     */
    public int getNVertices() {
        return vertices.length;
    }

    /**
     * Return (sub)profile's vertices. The real source and real receiver ground points are replaced by,
     * respectively, the source and the receiver
     *
     * @param iStart index of the first segment
     * @param iEnd index of the last segment
     * @return (sub)profile vertices
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
     * Return the coordinate of the ground vertex under the real source if i==0, the ground vertex under the real
     * receiver if i==nVertices-1 and vertices[i] otherwise
     *
     * @param i index of the (secondary) source
     * @return coordinate of the source ground point
     */
    private Coordinate getGroundVertex(int i){
        Coordinate groundPoint = (Coordinate) vertices[i].clone();
        if (i == 0){
            groundPoint.setY(groundPoint.getY() - sourceHeight);
        } else if (i == getNVertices()-1){
            groundPoint.setY(groundPoint.getY() - receiverHeight);
        }
        return groundPoint;
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
        double y1 = getGroundVertex(iSeg).getY();
        double y2 = getGroundVertex(iSeg+1).getY();
        double distance = getGroundVertex(iSeg).distance(getGroundVertex(iSeg+1));
        double xn = (y1 - y2) / distance;
        double yn = (x2 - x1) / distance;
        double d = (x0 - x1) * xn + (y0 - y1) * yn;
        double xi = x0 - 2 * d * xn ;
        double yi = y0 - 2 * d * yn;
        return new Coordinate(xi, yi);
    }

    /**
     * Return the abscissa of a point in a local coordinate system with the origin at the
     * normal projection of the source on the ground segment plane and the abscissa along the plane
     *
     * @param iPoint index of the point
     * @param iSeg index of the ground segment
     * @param iSource index of the (secondary) source
     * @return local abscissa (or local distance)
     */
    public double getLocalAbscissa(int iPoint, int iSeg, int iSource){
        double thetaPoint;
        thetaPoint = Angle.angleBetween(vertices[iPoint], getGroundVertex(iSeg + 1), getGroundVertex(iSeg));
        double thetaSource = Angle.angleBetween(vertices[iSource], getGroundVertex(iSeg+1), getGroundVertex(iSeg));
        return getGroundVertex(iSeg+1).distance(vertices[iSource]) * Math.cos(thetaSource)
                - getGroundVertex(iSeg+1).distance(vertices[iPoint]) * Math.cos(thetaPoint);
    }

    /**
     * Return the ordinate of a point in a local coordinate system with the origin at the
     * normal projection of the source on the ground segment plane and the abscissa along the plane
     * @param iPoint index of the point
     * @param iSeg index of the ground segment
     * @return local abscissa (or local height)
     */
    public double getLocalOrdinate(int iPoint, int iSeg) {
        double theta = Angle.angleBetweenOriented(vertices[iPoint], getGroundVertex(iSeg+1), getGroundVertex(iSeg));
        return getGroundVertex(iSeg+1).distance(vertices[iPoint]) * Math.sin(theta);
    }

    /**
     * Return the reflexion angle between source and receiver with respect to the normal on the ground
     * segment surface
     *
     * @param iSeg iSeg index of the ground segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return reflexion angle in radian
     */
    public double getReflexionAngle(int iSeg, int iSource, int iReceiver) throws IllegalArgumentException {
        if (isConvexSegment(iSeg, iSource, iReceiver)) {
            throw new IllegalArgumentException("Reflexion angle cannot be computed for convex segments");
        } else {

            Coordinate p = Intersection.intersection(getGroundVertex(iSeg), getGroundVertex(iSeg+1),
                    getImageVertex(iSource, iSeg), getVertex(iReceiver));
            Vector2D ps = new Vector2D(p, getVertex(iSource));
            Vector2D reflexionPlane = new Vector2D(getGroundVertex(iSeg), getGroundVertex(iSeg+1));
            return reflexionPlane.angle(ps) - Math.PI / 2;
        }
    }

    /**
     * Check if a segment is convex according to Harmonoise definition
     * Ref: Figure 3 and "Geometry" subsection of section 2.4.1 from Salomons et al.
     *
     * @param iSeg index of the ground segment
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return true if the segment is convex
     */
    public boolean isConvexSegment(int iSeg, int iSource, int iReceiver){
        double localSourceHeight = getLocalOrdinate(iSource, iSeg);
        double localReceiverHeight = getLocalOrdinate(iReceiver, iSeg);
        return localSourceHeight < 0 || localReceiverHeight < 0;
    }

    /**
     * Determine whether a ground (sub)profile contains convex segment or not
     * Ref: Figure 3 and "Geometry" subsection of section 2.4.1 from Salomons et al.
     *
     * @param iSource index of the (secondary) source
     * @param iReceiver index of the (secondary) receiver
     * @return true if the (sub)profile contains at least one convex segment
     */
    public boolean hasConvexSegment(int iSource, int iReceiver){
        boolean isConvex = false;
        // Loop on segments
        for (int i = iSource; i < iReceiver - 1; i++) {

            if (isConvexSegment(i, iSource, iReceiver)){
                isConvex = true;
                break;
            }
        }
        return isConvex;
    }

    /**
     * Get ground segment normalized ground impedance
     *
     * @param iSeg index of the ground segment
     * @param frequency frequency axis
     * @return ground impedance
     */
    public Complex[] getGroundImpedance(int iSeg, List<Double> frequency) {
        Complex[] groundImpedance = new Complex[frequency.size()];
        for (int i = 0; i < frequency.size(); i++) {
            groundImpedance[i] = SurfaceAbsorption.computeSurfaceImpedance(flowResistivity[iSeg], frequency.get(i));
        }
        return groundImpedance;
    }

    /**
     * Set the frequency dependant reflectionCoefficient for segment iSeg
     *
     * @param iSeg index of the ground segment
     * @param reflectionCoefficient spherical-wave reflection coefficient
     */
    public void setReflectionCoefficient(int iSeg, Complex[] reflectionCoefficient) {
        this.reflectionCoefficient[iSeg] = reflectionCoefficient;
    }

    /**
     * Return the reflection coefficient at a given frequency
     * @param iSeg index of the ground segment
     * @param iFreq index of the frequency
     * @return reflection coefficient
     */
    public Complex getReflectionCoefficient(int iSeg, int iFreq) {
        return reflectionCoefficient[iSeg][iFreq];
    }
}
