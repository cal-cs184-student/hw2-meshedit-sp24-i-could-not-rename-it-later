#include "student_code.h"
#include "CGL/vector3D.h"
#include "mutablePriorityQueue.h"
#include <complex>
#include <vector>

using namespace std;

namespace CGL
{

	/**
	 * Evaluates one step of the de Casteljau's algorithm using the given points and
	 * the scalar parameter t (class member).
	 *
	 * @param points A vector of points in 2D
	 * @return A vector containing intermediate points or the final interpolated vector
	 */
	std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
	{ 
		// TODO Part 1.
		std::vector<Vector2D> retval;
		for (int i = 0; i < points.size() - 1; i++) {
			retval.push_back(points[i] * (1 - t) + points[i + 1] * t);
		}

		return retval;
	}

	/**
	 * Evaluates one step of the de Casteljau's algorithm using the given points and
	 * the scalar parameter t (function parameter).
	 *
	 * @param points    A vector of points in 3D
	 * @param t         Scalar interpolation parameter
	 * @return A vector containing intermediate points or the final interpolated vector
	 */
	std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
	{
		// TODO Part 2.
		std::vector<Vector3D> retval;
		for (int i = 0; i < points.size() - 1; i++) {
			retval.push_back(points[i] * (1 - t) + points[i + 1] * t);
		}

		return retval;
	}

	/**
	 * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
	 *
	 * @param points    A vector of points in 3D
	 * @param t         Scalar interpolation parameter
	 * @return Final interpolated vector
	 */
	Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
	{
		// TODO Part 2.
		std::vector<Vector3D> tmp = points;

		while (tmp.size() > 1) {
			tmp = evaluateStep(tmp, t);
		}

		return tmp[0];
	}

	/**
	 * Evaluates the Bezier patch at parameter (u, v)
	 *
	 * @param u         Scalar interpolation parameter
	 * @param v         Scalar interpolation parameter (along the other axis)
	 * @return Final interpolated vector
	 */
	Vector3D BezierPatch::evaluate(double u, double v) const 
	{  
		// TODO Part 2.
		std::vector<Vector3D> tmp;
		for (int i = 0; i < controlPoints.size(); i++) {
			tmp.push_back(evaluate1D(controlPoints[i], u));
		}
		return evaluate1D(tmp, v);
	}


	Vector3D Vertex::normal( void ) const
	{
		// TODO Part 3.
		// Returns an approximate unit normal at this vertex, computed by
		// taking the area-weighted average of the normals of neighboring
		// triangles, then normalizing.
		
		std::vector<Vector3D> outer_vertices;
		std::vector<Vector3D> cross_products;
		Vector3D retval = Vector3D();

		HalfedgeCIter h = halfedge();         // get the outgoing half-edge of the vertex
		do {
			HalfedgeCIter h_twin = h->twin(); // get the opposite half-edge
			VertexCIter v = h_twin->vertex(); // vertex is the 'source' of the half-edge, so
											  // h->vertex() is v, whereas h_twin->vertex()
											  // is the neighboring vertex
			outer_vertices.push_back(v->position);

			h = h_twin->next();               // move to the next outgoing half-edge of the vertex
		} while(h != halfedge());             // keep going until we are back where we were
		
		for (int i = 0; i < degree(); i++) {
			Vector3D edge1 = outer_vertices[i % degree()] - position;
			Vector3D edge2 = outer_vertices[(i + 1) % degree()] - position;

			cross_products.push_back(cross(edge2, edge1));
		}

		for (int i = 0; i < degree(); i++) {
			retval += cross_products[i];
		}
		
		return retval.unit();
	}

	EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
	{
		// TODO Part 4.
		// This method should flip the given edge and return an iterator to the flipped edge.
		return EdgeIter();
	}

	VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
	{
		// TODO Part 5.
		// This method should split the given edge and return an iterator to the newly inserted vertex.
		// The halfedge of this vertex should point along the edge that was split, rather than the new edges.
		return VertexIter();
	}



	void MeshResampler::upsample( HalfedgeMesh& mesh )
	{
		// TODO Part 6.
		// This routine should increase the number of triangles in the mesh using Loop subdivision.
		// One possible solution is to break up the method as listed below.

		// 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
		// and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
		// a vertex of the original mesh.

		// 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.

		// 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
		// information about which subdivide edges come from splitting an edge in the original mesh, and which edges
		// are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
		// the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)

		// 4. Flip any new edge that connects an old and new vertex.

		// 5. Copy the new vertex positions into final Vertex::position.

	}
}
