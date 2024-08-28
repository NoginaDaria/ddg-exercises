// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    Eigen::SparseMatrix<size_t> A0(mesh->nEdges(), mesh->nVertices());

    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh->nEdges() * 2);

    size_t i = 0;

    for (Edge e : mesh->edges()) {
        i = geometry->edgeIndices[e];
        tripletList.push_back(T(i, geometry->vertexIndices[e.firstVertex()], 1));
        tripletList.push_back(T(i, geometry->vertexIndices[e.secondVertex()], 1));
    }

    A0.setFromTriplets(tripletList.begin(), tripletList.end());
    return A0;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    Eigen::SparseMatrix<size_t> A1(mesh->nFaces(), mesh->nEdges());

    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh->nFaces() * 3);

    size_t i = 0;

    for (Face f : mesh->faces()){
        i = geometry->faceIndices[f];
        tripletList.push_back(T(i, geometry->edgeIndices[f.halfedge().edge()], 1));
        tripletList.push_back(T(i, geometry->edgeIndices[f.halfedge().next().edge()], 1));
        tripletList.push_back(T(i, geometry->edgeIndices[f.halfedge().next().next().edge()], 1));
    }


    A1.setFromTriplets(tripletList.begin(), tripletList.end());
    return A1;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    Vector<size_t> vertices = Vector<size_t>::Zero(mesh->nVertices());
    for (size_t idx : subset.vertices) {
        vertices[idx] = 1;
    }
    return vertices;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    Vector<size_t> edges = Vector<size_t>::Zero(mesh->nEdges());
    for (size_t idx : subset.edges) {
        edges[idx] = 1;
    }
    return edges;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    Vector<size_t> faces = Vector<size_t>::Zero(mesh->nFaces());
    for (size_t idx : subset.faces) {
        faces[idx] = 1;
    }
    return faces;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    MeshSubset starSubset = subset.deepCopy();
    starSubset.addEdges(copyFromVector(this->A0 * this->buildVertexVector(starSubset)));
    starSubset.addFaces(copyFromVector(this->A1 * this->buildEdgeVector(starSubset)));
    return starSubset;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    MeshSubset closureSubset = subset.deepCopy();
    closureSubset.addEdges(copyFromVector(this->A1.transpose() * this->buildFaceVector(closureSubset)));
    closureSubset.addVertices(copyFromVector(this->A0.transpose() * this->buildEdgeVector(closureSubset)));
    return closureSubset;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    MeshSubset cl_st = closure(star(subset));
    cl_st.deleteSubset(star(closure(subset)));
    return cl_st;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    return subset.equals(closure(subset));
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    if (subset.faces.empty() && subset.edges.empty()) return 0;

    if (!subset.faces.empty()) {
        MeshSubset closured_subset({}, {}, subset.faces);
        return subset.equals(closure(closured_subset)) ? 2: -1;
    } else {
        MeshSubset closured_subset({}, subset.edges, {});
        return subset.equals(closure(closured_subset)) ? 1 : -1;
        
    }
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    MeshSubset boundary_subset;

    if (! isPureComplex(subset)) return boundary_subset;

    MeshSubset closureSubset = closure(subset);

    if (!subset.faces.empty()) {
        Vector<size_t> edgeVector = this->A1.transpose()*this->buildFaceVector(subset);
        for (size_t idx = 0; idx < (size_t)edgeVector.size(); idx++) {
            if (edgeVector[idx] == 1) {
                boundary_subset.addEdge(idx);
                
                Edge e = mesh->edge(idx);
                boundary_subset.addVertex(e.firstVertex().getIndex());
                boundary_subset.addVertex(e.secondVertex().getIndex());
            }
        }
    }

    if (!subset.edges.empty()) {
        Vector<size_t> vertexVector = this->A0.transpose()*this->buildEdgeVector(subset);
        for (size_t idx = 0; idx < (size_t)vertexVector.size(); idx++) {
            if (vertexVector[idx] == 1) {
                boundary_subset.addVertex(idx);
            }
        }
    }

    return boundary_subset;
}