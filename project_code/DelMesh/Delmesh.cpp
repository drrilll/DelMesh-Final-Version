//============================================================================
// Name        : Delmesh.cpp
// Author      : darryl
// Version     :
// Copyright   : sure
// Description : Delaunay mesh generator
//============================================================================

#include <iostream>
#include "Delmesh.h"
#include <functional>
#include <math.h>
#include <OpenMesh/Core/Mesh/PolyMeshT.hh>

using namespace std;



DelMesh::DelMesh(int ds, int st, string input, string output,  bool median){

    if (!OpenMesh::IO::read_mesh(mesh, input)){
        cout << "read error"<<endl;
    }


    this->output = output;
    score_type = st;
    // Add the samples property, which is a vector of potential vertex
    // sites along an edge.
    mesh.add_property(samples);

    //add indicator variable for non-Delaunay edges
    mesh.add_property(is_NDE);

    //add indicator for a flippable edge
    mesh.add_property(is_flippable);

    //requesting face colors, although currently not working
    mesh.request_face_colors();
    mesh.request_face_status();
    mesh.request_edge_status();
    mesh.request_vertex_status();

    //make the constants we need
    make_constants();

    cout<<"pv: "<<pv<<endl;
    cout<<"pe: "<<pe<<endl;

    //make non-Delaunay edges data structure
    q = new my_p_queue(&mesh, ds);

    //storing flippable NDE for processing
    flips = new my_p_queue(&mesh, 2);

    //a helper class for geometry functions with access to DelMesh state
    g2d = new Geom_2D(&mesh, &samples, &is_flippable, &is_NDE, st, median);


}

DelMesh::~DelMesh(){
    //write the mesh to file
    cout<<"writing "<<output<<endl;
    if (!OpenMesh::IO::write_mesh(mesh, output, writeOptions)){
        cout << "write error"<< endl;
    }
    delete q;
    delete g2d;
    delete flips;

    Mesh::EdgeIter eBegin, eEnd, eIt;
    eBegin = mesh.edges_begin();
    eEnd = mesh.edges_end();

    //delete all the sample points
    for (eIt = eBegin; eIt != eEnd; eIt++){
        if (!mesh.property(is_flippable, eIt) &&!mesh.is_boundary(eIt)){
            delete(mesh.property(samples, eIt));
        }
    }
}

/**
 * We can only call this once, otherwise we will have a memory leak.
 * @brief DelMesh::process_mesh
 */
void DelMesh::process_mesh(){
    // find all the NDE's
    cout<<"finding nde's"<<endl;

    // important that this is only called once, otherwise memory leak
    find_nd_edges();

    cout<<"done finding nde's"<<endl;

    cout<<"making the mesh"<<endl;
    make_Delaunay_mesh();
    cout<<"done making mesh, entering loop. Edges left to process:"<<endl;

    refind_nd_edges();

    while(q->size()>0){
        make_Delaunay_mesh();
        refind_nd_edges();
    }
    cout<<endl;

}

/**
 * Find all NDE's, but without making new sample points, and
 * store them in the appropriate data structure for processing
 * @brief DelMesh::refind_nd_edges
 */
void DelMesh::refind_nd_edges(){
    Mesh::EdgeIter eBegin, eEnd, eIt;
    eBegin = mesh.edges_begin();
    eEnd = mesh.edges_end();

    int count = 0;

    for (eIt = eBegin; eIt != eEnd; eIt++){
        //add sample points to the edge
        if (!mesh.is_boundary(*eIt)){
            if (is_nd_edge(*eIt, false)){
                if (mesh.property(is_flippable, *eIt)==FALSE){
                    if (mesh.property(samples, *eIt)->size()==0){
                        cout << "bad edge, length "<<mesh.calc_edge_length(*eIt)<<endl;
                        count ++;
                    }else{
                        q->push(*eIt);
                    }

                }else{
                    //cout<<"flipping"<<endl;
                    mesh.flip(*eIt);
                }
            }
        }
    }
    //cout<<"bad edges: "<<count<<" stack edges: "<<q->size()<<endl;
    cout<<q->size()<<" "<<flush;
}

/**
     * @brief make_Delaunay_mesh
     * All current NDE's have been found, now we process them and create the Delaunay mesh.
     * There is a lot of redundant code, that is, code that was copy and pasted. I thought of
     * refactoring into a function, but there would be a lot of flags to pass in for different
     * scenarios. For instance, newly created edges are processed twice, once for each half edge,
     * but we only want to transfer samples over once. Also, we don't want to add an edge multiple
     * times to a processing data structure. As such, the less complex approach seemed to be
     * copy and pasting the code and modifying it.
     */
void DelMesh::make_Delaunay_mesh(){

    Mesh::EdgeHandle eh;
    vector<Mesh::Point>* samps;
    int index;
    Mesh::Point point;
    Mesh::Point p1, p2;
    Mesh::HalfedgeHandle heh1, heh2;
    Mesh::VertexHandle to, from, vh1, vh2, mid, v;
    Mesh::FaceHandle fh1, fh2;
    vector<Mesh::VertexHandle> vec;

    int dummy_count = 0;
    int count = 0;
    while(!q->empty()){
        //cout<<"counting iterations..."<<endl;
        count ++;
        //if (count >1) return;
        //cout<<"edges: "<<mesh.n_edges()<<" count: "<<count<<" Stack size: "<<q->size()<<endl;
        eh = q->top();
        q->pop();

        //find the right sample point
        index = g2d->get_sample_point(eh);

        if (index == -1){
            //There is no provision for this, because theoretically it should
            //never happen. But it does.
            continue;
            //we pretend it doesn't. All the edges are processed in the end, so I
            //am not sure what happens here.
        }

        // remember to delete this when we are done
        samps = mesh.property(samples, eh);

        point = samps->at(index);

        //cout<<"next message samps index point: "<<samps->at(index)<<endl;

        heh1 = mesh.halfedge_handle(eh, 0);
        heh2 = mesh.halfedge_handle(eh, 1);

        to = mesh.to_vertex_handle(heh1);
        from = mesh.from_vertex_handle(heh1);

        vh1 = mesh.opposite_vh(heh1);
        vh2 = mesh.opposite_vh(heh2);


        //grab the two faces we wish to delete
        fh1 = mesh.face_handle(heh1);
        fh2 = mesh.face_handle(heh2);

        mesh.delete_face(fh1, false);
        mesh.delete_face(fh2, false);
        mesh.garbage_collection();

        /*
         * Add the new vertex and the new faces to the mesh, while simultaneously
         * marking flippable edges, splitting up sample points, and checking and
         * adding new edges
         */

        mid = mesh.add_vertex(point);

        /* From here we add 4 faces. However, we need access to the edges, and the
         * only way to get them is using an iterator, which does not start at
         * any particular edge. So we have to test each edge using, in this case,
         * vertex handles.
         */

        Mesh::EdgeHandle eh;

        /************* Face 1 *******************/
        vec.clear();
        vec.push_back(from);
        vec.push_back(mid);
        vec.push_back(vh1);
        Mesh::FaceHandle fh1 = mesh.add_face(vec);

        /************* Face 2 *******************/
        vec.clear();
        vec.push_back(from);
        vec.push_back(vh2);
        vec.push_back(mid);

        Mesh::FaceHandle fh2 = mesh.add_face(vec);

        /************* Face 3 *******************/
        vec.clear();
        vec.push_back(to);
        vec.push_back(vh1);
        vec.push_back(mid);

        Mesh::FaceHandle fh3 = mesh.add_face(vec);

        /************* Face 4 *******************/
        vec.clear();
        vec.push_back(mid);
        vec.push_back(vh2);
        vec.push_back(to);
        Mesh::FaceHandle fh4 = mesh.add_face(vec);

        //add in the samples in the split edge
        Mesh::FaceHalfedgeIter fhIt = mesh.fh_iter(fh1);

        Mesh::HalfedgeHandle heh;
        Mesh::VertexHandle v1, v2;



        //iterate over the edges, figure out what we got and how to handle it
        for (int i = 0; i < 3; i ++){
            eh = mesh.edge_handle(*fhIt);
            v = mesh.from_vertex_handle(*fhIt);

            // edge (from, mid)
            if (v == from){
                // this is a new edge made from the old edge, we need new samples
                mesh.property(samples, eh) = new vector<Mesh::Point>();
                //we want to put the samples on, but we have to
                //put the correct ones
                if (is_on_edge(samps->at(0), eh)){
                    for (int i = 0; i < index; i++){
                        mesh.property(samples, eh)->push_back(samps->at(i));
                    }
                }else if (is_on_edge(samps->at(samps->size()-1), eh)){
                    for (int i = index +1; i < samps->size(); i++){
                        mesh.property(samples, eh)->push_back(samps->at(i));
                    }
                }else{
                    dummy_count ++;
                }

                // edge (mid, vh1)
            }else if (v == mid){
                mesh.property(is_flippable, eh) = TRUE;
                // edge (vh1, from)
            }else if (v == vh1){
                if (mesh.property(is_NDE, eh)==1){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        flip_later(eh);
                    }
                }else if (is_nd_edge(eh, false) && !mesh.is_boundary(eh)){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        flip_later(eh);
                    }else{
                        mesh.property(is_NDE, eh) = TRUE;
                        //cout<<"pushing vh1, from"<<endl;
                        //g2d->output_point(eh);
                        q->push(eh);
                    }
                }

            }else{
                cout<<"damn"<<endl;
            }
            fhIt++;
        }


        //Face 2
        fhIt = mesh.fh_iter(fh2);

        //iterate over the edges, figure out what we got and how to handle it
        //we have already added samples above, so we skip that step
        for (int i = 0; i < 3; i ++){
            eh = mesh.edge_handle(*fhIt);
            v = mesh.from_vertex_handle(*fhIt);
            if (v==mid){
                //nothing to be done
            }

            // edge (mid, from)
            if (v == vh2){
                mesh.property(is_flippable, eh) = TRUE;
                // edge (from, vh2)
            }else if (v == from){
                if (mesh.property(is_NDE, eh)==TRUE){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        flip_later(eh);
                    }
                }else if (is_nd_edge(eh, false) && !mesh.is_boundary(eh)){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        flip_later(eh);
                    }else{
                        mesh.property(is_NDE, eh) = TRUE;
                        q->push(eh);
                    }
                }

            }
            fhIt++;
        }


        //Face 3
        fhIt = mesh.fh_iter(fh3);

        //iterate over the edges, figure out what we got and how to handle it
        for (int i = 0; i < 3; i ++){
            eh = mesh.edge_handle(*fhIt);
            v = mesh.from_vertex_handle(*fhIt);

            // edge (mid, to)
            if (v == mid){
                mesh.property(samples, eh) = new vector<Mesh::Point>();
                //we want to put the samples on, but we have to
                //put the correct ones
                if (is_on_edge(samps->at(0), eh)){
                    for (int i = 0; i < index; i++){
                        mesh.property(samples, eh)->push_back(samps->at(i));
                    }
                }else if (is_on_edge(samps->at(samps->size()-1), eh)){
                    for (int i = index +1; i < samps->size(); i++){
                        mesh.property(samples, eh)->push_back(samps->at(i));
                    }
                }else{
                    dummy_count ++;
                }

                // edge (vh1, mid) has been handled
                // edge (to, vh1)
            }else if (v == to){
                if (mesh.property(is_NDE, eh)==TRUE){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        flip_later(eh);
                    }
                }else if (is_nd_edge(eh, false) && !mesh.is_boundary(eh)){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        flip_later(eh);
                    }else{
                        mesh.property(is_NDE, eh) = TRUE;
                        q->push(eh);
                    }
                }

            }
            fhIt++;
        }


        //Face 4
        fhIt = mesh.fh_iter(fh4);

        //iterate over the edges, figure out what we got and how to handle it
        for (int i = 0; i < 3; i ++){
            eh = mesh.edge_handle(*fhIt);
            v = mesh.from_vertex_handle(*fhIt);

            if (v==to){
                //nothing to be done
            }

            // edge (to, mid) has been handled
            // edge (mid, vh2) has been handled
            // edge (vh2, to)
            if (v == vh2){
                if (mesh.property(is_NDE, eh)==TRUE){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        flip_later(eh);
                    }
                }else if (is_nd_edge(eh, false) && !mesh.is_boundary(eh)){
                    if (mesh.property(is_flippable, eh)==TRUE){
                        flip_later(eh);
                    }else{
                        mesh.property(is_NDE, eh) = TRUE;
                        q->push(eh);
                    }
                }

            }

            fhIt++;
        }

        //cout<<"stack size: "<<q->size()<<endl;

        //the circulators are done, so now flip the flippable edges
        while(flips->size()>0){
            mesh.flip(flips->top());
            flips->pop();
        }
        delete(samps);

    }
    //cout<<"dummy nodes: "<<dummy_count<<endl;

}

/**
 * @brief DelMesh::test_flip
 * @param eh
 *
 * Put the edges to be flipped in a stack to be flipped
 * after the circulators are done
 */
void DelMesh::flip_later(Mesh::EdgeHandle eh){
    flips->push(eh);
}

/**
 * We cheat a bit here. We find the midpoint of the line. Then see if the
 * vertex is within half the length of the midpoint.
 * @brief DelMesh::is_on_line
 * @param v
 * @param eh
 * @return
 */
bool DelMesh::is_on_edge(Mesh::Point &v, Mesh::EdgeHandle &eh)
{
    Mesh::Scalar length = mesh.calc_edge_length(eh)/2;
    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
    Mesh::Point to, from;
    to = mesh.point(mesh.to_vertex_handle(heh));
    from = mesh.point(mesh.from_vertex_handle(heh));
    if ((equals(v,to))||(equals(v,from))){
        return false;
    }

    to -= from;
    to.normalize();
    to *= length;
    to += from;
    bool answer = (g2d->distance3d(to, v)<= length);
    return answer;
}

bool DelMesh::equals(Mesh::Point p1, Mesh::Point p2){
    bool answer = ((p1[0]==p2[0])&&(p1[1]==p2[1])&&(p1[2]==p2[2]));
    return answer;
}


/**
 * @brief DelMesh::make_constants
 *
 * We make the constants that are used to determine the number and spacing of the sample
 * points. They are based on minimum angle and ratio minimum:maximum edge lengths
 */
void DelMesh::make_constants(){
    /*
         * Find the shortest edge, and the minimum angle
         */
    Mesh::VertexIter vIt, vBegin, vEnd;
    Mesh::VertexOHalfedgeIter veIt, veBegin, veEnd;;

    Mesh::Scalar min_angle = 7;
    Mesh::Scalar min_edge = 10000000;

    vBegin = mesh.vertices_begin();
    vEnd = mesh.vertices_end();
    int count = 0, c2 = 0;
    Mesh::Scalar length = 0, angle,maxl = 0;
    for (vIt = vBegin; vIt != vEnd; ++vIt){
        veBegin = mesh.voh_begin(*vIt);
        veEnd = mesh.voh_end(*vIt);

        for (veIt = veBegin; veIt != veEnd; ++veIt){
            c2++;
            length = mesh.calc_edge_length(*veIt);
            if (length > maxl){maxl = length;}
            if (length<min_edge && length != 0){
                min_edge = length;
            }
            angle = mesh.calc_sector_angle(*veIt);
            if (angle < 0){angle = - angle;}
            if (angle < min_angle&&angle != 0){
                min_angle = angle;
            }
            if (angle ==0){
                count++;
            }
        }
    }
    cout<<"zero angles: "<< count <<endl;
    cout<<"total angles: "<< c2 <<endl;
    cout<<"min angle: "<< min_angle <<endl;
    cout << "max length: "<< maxl<<endl;
    cout << "min length " << min_edge<<endl;


    pv = min_edge*sin(min_angle)/(0.5 + sin(min_angle));
    pe = min_edge/2;

    if (pe<pv){pv = pe;}

    pe = 2 * pv * sin(min_angle);

    double k = maxl/(min_edge*sin(min_angle)*sin(min_angle));

    cout <<"K: "<<k<<endl;

    double kk = maxl/ pe;

    cout <<"KK: "<<kk<<endl;
}


/*
     * We have been given a non-Delaunay edge, and now we will
     * add sample points to it.
     */
void DelMesh::make_sample_points(Mesh::EdgeHandle ehandle){
    /*
         * Get the two vertices associated with the edge, calculate and add samples
         * in a specified order. We always add the first sample point; after that
         * we check. For a definition of sample point, see the associated paper.
         */
    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(ehandle, 1);
    Mesh::Point to, from;
    //cout<<"getting to and from vertices"<<endl;
    to = mesh.point(mesh.to_vertex_handle(hedge));
    from = mesh.point(mesh.from_vertex_handle(hedge));

    double length = mesh.calc_edge_length(hedge);
    Mesh::Normal unit(mesh.calc_edge_vector(hedge));
    unit.normalize();

    double total = pv;
    Mesh::Point samp(unit);
    samp *= total;
    samp += from;

    mesh.property(samples, ehandle)->push_back(samp);

    // We have our first sample in. If this is the shortest edge, it could
    // be the only sample (highly unlikely though).
    if (length <= 2 * pv){
        return;
    }

    // Otherwise keep adding samples
    double mark = length - (pv + pe);
    while (total < mark){
        total += pe;
        samp = (unit * total) + from;
        mesh.property(samples, ehandle)->push_back(samp);
    }

    // We add one more sample at length *pv from the "to" vertex
    samp = (unit * (length - pv))+from;
    mesh.property(samples, ehandle)->push_back(samp);
    // all samples have been added

}


/*
     * Run some tests to make sure the samples generated are accurate
     */
void DelMesh::test_samples(){
    Mesh::EdgeIter edge = mesh.edges_begin();
    make_sample_points(*edge);
    vector<Mesh::Point>* samps;
    cout<<"getting test samples"<<endl;
    samps = (mesh.property(samples, *edge));
    cout<<"got test samples"<<endl;
    double length = mesh.calc_edge_length(*edge);
    int num_samples = 2;
    length -= 2*pv;

    length /= pe;
    num_samples += floor(length);
    cout << "expected samples: "<<num_samples<<endl;
    cout << "samples generated: "<<samps->size()<<endl;

    //now we want to check that they are all on the same line.
    //They should all output the same unit vector as the edge itself
    Mesh::Normal unit(mesh.calc_edge_vector(*edge));
    unit.normalize();
    cout<<endl<<"Normal vector: "<<unit<<endl;

    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(*edge, 1);

    Mesh::Point from = mesh.point(mesh.from_vertex_handle(hedge));
    Mesh::Point point(0,0,0);
    vector<Mesh::Point>::iterator it, it2;
    for (it = samps->begin(); it < samps->end(); it++){
        point = *it - from;
        point.normalize();
        cout<<endl<<"Normal vector: "<<unit<<endl;
        cout<<"Sample vector: "<<point<<endl;
        cout<<"Before normal: "<<*it<<endl;
    }

    edge ++; edge ++; edge ++;
    make_sample_points(*edge);
    vector<Mesh::Point>* sam;
    sam = (mesh.property(samples, *edge));
    it = samps->begin(); it2 = sam->begin();
    for (int i = 0; i < 10; i++){
        cout << "Point 1: "<<*it++<<endl;
        cout << "Point 2: "<<*it2++<<endl;

    }


}



/**
     * We check the angles opposite the edge to see if they
     * sum to > pi.
     */
bool DelMesh::is_nd_edge(Mesh::EdgeHandle edge, bool output = false){

    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(edge, 1);
    Mesh::HalfedgeLoopIter hIt = mesh.hl_begin(hedge);

    //we want the angle opposite this half-edge, which is the
    //sector angle of the next half-edge
    hIt++;
    Mesh::Scalar angle1 = mesh.calc_sector_angle(*hIt);

    // we need absolute values
    if (angle1 <0){
        angle1 = -angle1;
    }

    //repeat for the opposite half-edge to get the opposite angle
    Mesh::HalfedgeHandle he_opposite_handle = mesh.opposite_halfedge_handle(hedge);
    hIt = mesh.hl_begin(he_opposite_handle);
    hIt++;
    Mesh::Scalar angle2 = mesh.calc_sector_angle(*hIt);
    if (angle2 <0){
        angle2 = -angle2;
    }
    if (output&&(angle1+angle2)>M_PI){
        cout<<"angle: "<<angle1+angle2<<endl;
    }
    return ((angle1+angle2)>M_PI);

}
/*
     * Find all non-Delaunay edges and put them in a data structure for processing
     */
void DelMesh::find_nd_edges(){
    Mesh::EdgeIter eBegin, eEnd, eIt;
    eBegin = mesh.edges_begin();
    eEnd = mesh.edges_end();

    int count = 0;

    for (eIt = eBegin; eIt != eEnd; eIt++){
        count ++;
        //and now we are managing memory
        mesh.property(samples, *eIt) = new vector<Mesh::Point>();
        //add sample points to the edge
        if (!mesh.is_boundary(*eIt)){
            make_sample_points(*eIt);
            if (is_nd_edge(*eIt, false)){
                //indicate that the edge is non-Delaunay
                mesh.property(is_NDE, *eIt) = 1;
                //add the edge to our data structure of current NDE's
                q->push(*eIt);
            }
        }
    }
    cout<<"# NDE's: "<<q->size()<<endl;
    cout << "# edges: "<< count << endl;

}

Mesh* DelMesh::getMesh(){
    return &mesh;
}

Prop* DelMesh::getSamples(){
    return &samples;
}


/*
     * Test the priority queue. Presumably the NDE's have been added. Now
     * pop them off and output the lengths
     */
void DelMesh::test_pq(){
    Mesh::EdgeIter eBegin, eEnd, eIt;
    eBegin = mesh.edges_begin();
    eEnd = mesh.edges_end();

    my_p_queue stack(&mesh, 2);
    my_p_queue queue(&mesh, 1);

    /*
         * test it on 100 edges
         */

    eIt = eBegin;
    Mesh::Scalar length = 0;
    for (int i = 0; i < 100; i++, eIt++){
        q->push(*eIt);
        queue.push(*eIt);
    }

    for (int i = 0; i < 100; i++, eIt++){
        length = mesh.calc_edge_length(q->top());
        stack.push(q->top());
        queue.push(q->top());
        q->pop();
        cout << "length "<<i<<": "<<length<<endl;
    }

    for (int i = 0; i < 100; i++, eIt++){
        length = mesh.calc_edge_length(stack.top());
        stack.pop();
        cout << "stack length "<<i<<": "<<length<<endl;
    }

    for (int i = 0; i < 100; i++, eIt++){
        length = mesh.calc_edge_length(queue.top());
        stack.push(queue.top());
        queue.pop();
        cout << "queue length "<<i<<": "<<length<<endl;
    }

    for (int i = 0; i < 100; i++, eIt++){
        length = mesh.calc_edge_length(stack.top());
        stack.pop();
        cout << "stack length "<<i<<": "<<length<<endl;
    }


}

/*
     * Testing the face colors. Spoiler alert, they don't work. That is when
     * I try to load the file into a viewer after it throws and error.
     */
void DelMesh::test_face_colors(){
    if (!OpenMesh::IO::read_mesh(mesh, "ateneav.obj")){
        cout << "read error"<<endl;
    }
    if (!OpenMesh::IO::write_mesh(mesh, "ateneav2.obj", writeOptions)){
        cout << "write error"<< endl;
    }

}


int DelMesh::test_2D_flattening(){
    Mesh mesh;
    // generate vertices
    Mesh::VertexHandle vhandle[8];
    vhandle[0] = mesh.add_vertex(Mesh::Point(-1, -1,  1));
    vhandle[1] = mesh.add_vertex(Mesh::Point( 1, -1,  1));
    vhandle[2] = mesh.add_vertex(Mesh::Point( 1,  1,  1));
    vhandle[3] = mesh.add_vertex(Mesh::Point(-1,  1,  1));
    vhandle[4] = mesh.add_vertex(Mesh::Point(-1, -1, -1));
    vhandle[5] = mesh.add_vertex(Mesh::Point( 1, -1, -1));
    vhandle[6] = mesh.add_vertex(Mesh::Point( 1,  1, -1));
    vhandle[7] = mesh.add_vertex(Mesh::Point(-1,  1, -1));
    // generate (triangle) faces
    std::vector<Mesh::VertexHandle>  face_vhandles;
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[2]);
    mesh.add_face(face_vhandles);


    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[3]);
    mesh.add_face(face_vhandles);

    //*********************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[5]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[4]);
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[5]);
    mesh.add_face(face_vhandles);

    /**************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[4]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[4]);
    face_vhandles.push_back(vhandle[5]);
    mesh.add_face(face_vhandles);

    /*******************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[6]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[2]);
    mesh.add_face(face_vhandles);

    /****************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[3]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[3]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[7]);
    mesh.add_face(face_vhandles);

    /**************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[3]);
    face_vhandles.push_back(vhandle[7]);
    mesh.add_face(face_vhandles);

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[4]);
    mesh.add_face(face_vhandles);

    // write mesh to cube.obj

    try
    {
        if ( !OpenMesh::IO::write_mesh(mesh, "cube.obj") )
        {
            std::cerr << "Cannot write mesh to file 'cube.obj'" << std::endl;
        }
    }
    catch( std::exception& x )
    {
        std::cerr << x.what() << std::endl;
    }


    Geom_2D g(&mesh, &samples, &is_flippable, &is_NDE, score_type, true);

    Mesh::EdgeIter eIt = mesh.edges_begin();

    Mesh::HalfedgeHandle hedge = mesh.halfedge_handle(*eIt, 1);
    Mesh::Point to, from;
    to = mesh.point(mesh.to_vertex_handle(hedge));
    from = mesh.point(mesh.from_vertex_handle(hedge));

    cout<<"point 1 x:"<<to[0]<<" y: "<<to[1]<<" z: "<<to[2]<<endl;
    cout<<"point 2 x:"<<from[0]<<" y: "<<from[1]<<" z: "<<from[2]<<endl;


    vector<Point_2D> p(g.test_flattening(*eIt));


    Mesh mesh2;
    // generate vertices
    for (int i = 0; i < 8; i++){
        vhandle[i] = mesh2.add_vertex( (p[i]));
        std::cout<<"point "<<i<<": "<<(p[i])<<endl;
    }

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[2]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[5]);
    face_vhandles.push_back(vhandle[1]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    /*********************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[7]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[5]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[6]);
    face_vhandles.push_back(vhandle[5]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    /**************/
    face_vhandles.clear();
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[1]);
    face_vhandles.push_back(vhandle[3]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    face_vhandles.clear();
    face_vhandles.push_back(vhandle[0]);
    face_vhandles.push_back(vhandle[2]);
    face_vhandles.push_back(vhandle[4]);
    mesh2.add_face(face_vhandles);
    cout<<"complex?"<<endl;

    // write mesh to flat_cube.obj
    try
    {
        if ( !OpenMesh::IO::write_mesh(mesh2, "flat_cube.obj") )
        {
            std::cerr << "Cannot write mesh to file 'cube.obj'" << std::endl;
            return 1;
        }
    }
    catch( std::exception& x )
    {
        std::cerr << x.what() << std::endl;
        return 1;
    }
    return 0;
}








