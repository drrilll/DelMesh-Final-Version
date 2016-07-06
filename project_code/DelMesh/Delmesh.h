/*
 * Delmesh.h
 *
 *  Created on: Mar 14, 2016
 *      Author: darryl
 */

#ifndef DELMESH_H_
#define DELMESH_H_

#include "priorityqueue.h"
#include "geom_2d.h"

class DelMesh
{

private:

    Mesh mesh;
    Prop samples;
    Delaunay_indicator is_NDE;
    Flippable is_flippable;
    Mesh::Scalar pv, pe;
    my_p_queue* q, *flips;
    OpenMesh::IO::Options writeOptions;
    Geom_2D* g2d;

    /*
     * Originally preprocessor values, but recompiling was
     * taking forever.
     */

    string output;
    int score_type;
    int TRUE = 1;
    int FALSE = 0;

public:

    DelMesh(int ds, int st, string input, string output, bool median);
    ~DelMesh();

    /**
     * @brief process_Mesh
     * Do everything.
     */
    void process_mesh();


    /**
     * @brief make_Delaunay_mesh
     * All NDE's have been found, now we process them and create the Delaunay mesh.
     */
    void make_Delaunay_mesh();

    void refind_nd_edges();

    void flip_later(Mesh::EdgeHandle eh);

    /*
     * make the constants pe and pv - see paper
     */
    void make_constants();

    //Mesh::Scalar get_intervals(Mesh::HalfedgeHandle heh, Mesh::Scalar &i1, Mesh::Scalar &i2);


    /*
     * We have been given a non-Delaunay edge, and now we will
     * add sample points to it.
     */
    void make_sample_points(Mesh::EdgeHandle ehandle);

    vector<Mesh::Point>* get_samples(Mesh::Edge edge);

    /*
     * Run some tests to make sure the samples generated are accurate
     */
    void test_samples();
    /**
     * We check the angles opposite the edge to see if they
     * sum to <= pi.
     */
    bool is_nd_edge(Mesh::EdgeHandle edge, bool output);

    void process_face(Mesh::VertexHandle &from_split_edge, Mesh::VertexHandle &from_new_edge, Mesh::VertexHandle &from_old_edge,
                      Mesh::FaceHalfedgeIter &fhIt, bool sample, int index);

    /*
     * Find all non-Delaunay edges and put them in a priority queue for
     * max length
     */
    void find_nd_edges();

    bool equals(Mesh::Point p1, Mesh::Point p2);

    bool is_on_edge(Mesh::Point &v, Mesh::EdgeHandle &eh);
    void test_face_colors();

    Mesh* getMesh();

    Prop* getSamples();

    /*
     * Test the priority queue
     */
    void test_pq();
    int test_2D_flattening();




};



#endif /* DELMESH_H_ */
