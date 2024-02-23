#include "student_code.h"
#include "CGL/vector3D.h"
#include "halfEdgeMesh.h"
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
		HalfedgeIter h0 = e0->halfedge();
		HalfedgeIter h1 = h0->twin();

		// Don't flip edges between boundaries.
		if (h0->face()->isBoundary() || h1->face()->isBoundary())
			return e0;

		// Find four edges that need to be changed.
		HalfedgeIter ad = h0->next();
		HalfedgeIter cb = h1->next();
		HalfedgeIter dc = ad->next();
		HalfedgeIter ba = cb->next();

		// Set half-edges of flipped edge appropariately.
		h0->setNeighbors(ba,
						 h1,
						 dc->vertex(),
						 h0->edge(),
						 h0->face());

		h1->setNeighbors(dc,
						 h0,
						 ba->vertex(),
						 h1->edge(),
						 h1->face());
		// Set other affected edges appropriately.
		dc->setNeighbors(cb, dc->twin(), dc->vertex(), dc->edge(), h1->face());
		ba->setNeighbors(ad, ba->twin(), ba->vertex(), ba->edge(), h0->face());
		ad->setNeighbors(h0, ad->twin(), ad->vertex(), ad->edge(), h0->face());
		cb->setNeighbors(h1, cb->twin(), cb->vertex(), cb->edge(), h1->face());
		
		// Set vertices appropriately.
		cb->vertex()->halfedge() = cb;
		ad->vertex()->halfedge() = ad;
		// Set faces appropriately.
		cb->face()->halfedge() = cb;
		ad->face()->halfedge() = ad;

		//printf("halfedge: %p\nnext: %p, twin: %p\nvertex: %p, edge: %p, face: %p\n\n",
		// 		h0, h0->next(), h0->twin(), h0->vertex(), h0->edge(), h0->face());
		/*check_for(h0);
		printf("\n");
		check_for(h0->next());
		printf("\n");
		check_for(h0->next()->next());
		printf("~~~~~~~~~~~~~~~~~~\n");
		check_for(h1);
		printf("\n");
		check_for(h1->next());
		printf("\n");
		check_for(h1->next()->next());
		printf("------------------\n");*/
		return e0;
	}

	VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
	{
		// TODO Part 5.
		// This method should split the given edge and return an iterator to the newly inserted vertex.
		// The halfedge of this vertex should point along the edge that was split, rather than the new edges.
		//
		HalfedgeIter h0 = e0->halfedge();
		HalfedgeIter h1 = e0->halfedge()->twin();

		VertexIter a = h1->vertex();
		VertexIter b = h1->next()->twin()->vertex();
		VertexIter c = h0->vertex();
		VertexIter d = h0->next()->twin()->vertex();

		FaceIter f0 = h0->face();
		FaceIter f1 = h1->face();

		// Get counterclockwise halfedges
		HalfedgeIter ad = h0->next(); HalfedgeIter dc = ad->next();
		HalfedgeIter cb = h1->next(); HalfedgeIter ba = cb->next();
		
		// Create m and assign its fields
		VertexIter m = newVertex();
		m->getVertex()->position = a->position / 2 + c->position / 2;
		m->halfedge() = h0;

		// Create new halfedges
		HalfedgeIter h2 = newHalfedge(); HalfedgeIter h3 = newHalfedge();
		HalfedgeIter h4 = newHalfedge(); HalfedgeIter h5 = newHalfedge();
		HalfedgeIter h6 = newHalfedge(); HalfedgeIter h7 = newHalfedge();

		// Handle halfedge twins
		h0->twin() = h1; h1->twin() = h0;
		h2->twin() = h3; h3->twin() = h2;
		h4->twin() = h5; h5->twin() = h4;
		h6->twin() = h7; h7->twin() = h6;

		// Assign halfedge nexts
		h7->next() = h0; h0->next() = ad; ad->next() = h7; 
		h1->next() = h2; h2->next() = ba; ba->next() = h1;
		h3->next() = h4; h4->next() = cb; cb->next() = h3;
		h5->next() = h6; h6->next() = dc; dc->next() = h5;

		// Assign even halfedge vertices  	// Assign odd halfedge vertices
		h0 -> vertex() = m;               	h1 -> vertex() = a; a -> halfedge() = h1;
		h2 -> vertex() = m;               	h3 -> vertex() = b; b -> halfedge() = h3;
		h4 -> vertex() = m;               	h5 -> vertex() = c; c -> halfedge() = h5;
		h6 -> vertex() = m;               	h7 -> vertex() = d; d -> halfedge() = h7;

		// Create new edges
		EdgeIter e1 = newEdge();
		EdgeIter e2 = newEdge();
		EdgeIter e3 = newEdge();

		// Assign each edge an halfedge     // Assign each halfedge an edge.
		e0->halfedge() = h0;                h0->edge() = e0; h1->edge() = e0;
		e1->halfedge() = h2;                h2->edge() = e1; h3->edge() = e1;
		e2->halfedge() = h4;                h4->edge() = e2; h5->edge() = e2;
		e3->halfedge() = h6;                h6->edge() = e3; h7->edge() = e3;

		// Create new faces
		FaceIter f2 = newFace();
		FaceIter f3 = newFace();

		// Assign faces to halfedges
		h7->face() = f0; h0->face() = f0; ad->face() = f0;
		h1->face() = f1; h2->face() = f1; ba->face() = f1;
		h3->face() = f2; h4->face() = f2; cb->face() = f2;
		h5->face() = f3; h6->face() = f3; dc->face() = f3;

		// Assign halfedges to faces
		f0->halfedge() = h0; f1->halfedge() = h2; f2->halfedge() = h4; f3->halfedge() = h6;

		return m;
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
