#include "geom_2d.h"
#include "triangle_properties.hpp"

Geom_2D::Geom_2D(Mesh *mesh, Prop *samples, Flippable *is_flippable, Delaunay_indicator* is_NDE, int score_type, bool median)
{
    this->mesh = mesh;
    this->samples = samples;
    this->is_flippable = is_flippable;
    this->is_NDE = is_NDE;
    this->score_type = score_type;
    this->median = median;
}

/**
 * The goal is to return the best sample point that we can split this
 * edge on. Returns the index of said sample point.
 */
int Geom_2D::get_sample_point(Mesh::EdgeHandle &ehandle){
    /*
     * The goal is to convert these coordinates into 2 dimensions
     * and do the necessary calculations there. The position of the
     * translated points is arbitrary, we shall merely calculate a
     * distance then use a unit vector to find the relevant sample point
     */

    vector<int> stack;

    Mesh::HalfedgeHandle heh = mesh->halfedge_handle(ehandle, 0);

    //starting point. We position our edge along the x-axis.

    vector<Point_2D> p;
    p.push_back(Point_2D(0.0, 0.0, 0.0));  //0

    p.push_back(Point_2D(mesh->calc_edge_length(heh), 0.0, 0.0)); //1

    p.push_back(get_2D_point(heh, p[0], p[1])); //2

    Mesh::HalfedgeLoopIter heIt;
    heIt = mesh->hl_begin(heh);

    heIt ++;
    p.push_back(get_2D_point(mesh->opposite_halfedge_handle(*heIt), p[2], p[1]));  //3

    heIt ++;
    p.push_back(get_2D_point(mesh->opposite_halfedge_handle(*heIt), p[0], p[2])); //4

    Mesh::HalfedgeHandle heh2 = mesh->opposite_halfedge_handle(heh);

    p.push_back(get_2D_point(heh2, p[1], p[0])); //5
    heIt = mesh->hl_begin(heh2);

    heIt ++;
    p.push_back(get_2D_point(mesh->opposite_halfedge_handle(*heIt), p[5], p[0])); //6

    heIt ++;
    p.push_back(get_2D_point(mesh->opposite_halfedge_handle(*heIt), p[1], p[5])); //7

    /*
     * We aim for simplicity and accuracy over speed in this section.
     * cc and r contain the circumcenter and the radius of the circumcircles
     * of the provided triangles.
     */
    double cc[2*6], r[6];

    get_circumcircle(p[0], p[2], p[5], &cc[0], r[0]);
    get_circumcircle(p[1], p[2], p[5], &cc[2], r[1]);
    get_circumcircle(p[1], p[2], p[3], &cc[4], r[2]);
    get_circumcircle(p[1], p[7], p[5], &cc[6], r[3]);
    get_circumcircle(p[0], p[2], p[4], &cc[8], r[4]);
    get_circumcircle(p[0], p[6], p[5], &cc[10], r[5]);

    vector<Mesh::Point>* samps;
    samps = (mesh->property(*samples, ehandle));
    double score;
    Point_2D p2d(0.0,0.0,0.0);
    int maxScore = 0;
    int maxIndex = 0;

    if (samps->size()==0){
        //this should never happen....but it does.
        return -1;

    }

    for (int i = 0; i < samps->size(); i++){
        //first two circles if the sample points are in both we
        //score 5
        //void Geom_2D::mesh_to_plane(Mesh::HalfedgeHandle heh, Mesh::Point &samp, Point_2D &p0, Point_2D &p1, Point_2D &dest){
        score = 0;
        mesh_to_plane(heh, (*samps)[i], (p[0]), (p[1]), p2d);

        if (score_type == 1){
            //cout<<"scoring"<<endl;
            if (distance2d(cc, p2d)<r[0]){
                score += 1.5;
                //cout <<score<<endl;
            }
            if (distance2d(cc+2, p2d)<r[1]){
                score += 1.5;
                //cout <<score<<endl;
            }
            //cout <<score<<endl;
        }else if (score_type == 0){
            if (distance2d(cc, p2d)<r[0]){
                score ++;
            }
            if (distance2d(cc+2, p2d)<r[1]){
                score ++;
            }
        }else if (score_type == 2){
            if ((distance2d(cc, p2d)<r[0])&&(distance2d(cc+2, p2d)<r[1])){
                score += 5;
            }
        }

        for (int j = 2; j < 6; j++){
            if (distance2d(cc+(2*j), p2d)>r[j]){
                score++;
            }
        }
        if (score>maxScore){
            maxScore = score;
            maxIndex = i;
            stack.clear();
            stack.push_back(i);
        }
        if (score == maxScore){
            stack.push_back(i);
        }

    }

    //this should almost never happen either. Of course sometimes ...
    if(samps->size()==0){
        cout<<"*************************"<<endl;
        cout<<"******************index: "<<maxIndex<<endl;
        cout<<"*************************"<<endl;
        cout<<"samples size: "<<samps->size()<<endl;
        cout<<" is flippable "<< mesh->property(*is_flippable, ehandle)<<endl;
        cout<<" is Delaunay "<< mesh->property(*is_NDE, ehandle)<<endl;
        output_point(ehandle);
        cout<<"length: "<<mesh->calc_edge_length(ehandle)<<endl;
        cout<<"max score "<<maxScore<<endl;

    }
    if(!median){
        return maxIndex;
    }

    int index = stack.size()/2;
    if (stack.size() == 0){
        //this is a problem, but I can't remember the last time it happened
        cout<<"***********************************"<<endl;
        cout<<"***********************************"<<endl;
        cout<<"***********************************"<<endl;
        cout<<"***********************************"<<endl;
        cout<<"***********************************"<<endl;

        return 0;
    }
    return stack.at(index);
}

void Geom_2D::output_point(Mesh::EdgeHandle eh){
    Mesh::Point p1 = mesh->point(mesh->to_vertex_handle(mesh->halfedge_handle(eh,0)));
    Mesh::Point p2 = mesh->point(mesh->from_vertex_handle(mesh->halfedge_handle(eh,0)));

    cout <<"point1: "<<p1<<endl;
    cout <<"point2: "<<p2<<endl;
}


/**
 * The assumption is that p0 corresponds to the from vertex. From there we
 * map the sample point into the 2D plane
 * @brief Geom_2D::mesh_to_plane
 * @param heh
 * @param samp
 * @param p0
 * @param p1
 */
void Geom_2D::mesh_to_plane(Mesh::HalfedgeHandle heh, Mesh::Point &samp, Point_2D &p0, Point_2D &p1, Point_2D &dest){

    Mesh::Point from = mesh->point(mesh->from_vertex_handle(heh));
    Mesh::Scalar dist = distance3d(from, samp);
    if (p0[0]<p1[0]){
        dest[0] = p0[0]+dist;
    }else{
        dest[0] = p0[0]-dist;
    }
    dest[1] = 0.0;
}

Mesh::Scalar Geom_2D::distance2d(double* p1, Point_2D p2){
    return sqrt(square(p1[0]-p2[0]) + square(p1[1]-p2[1]));
}

Mesh::Scalar Geom_2D::distance3d(Mesh::Point p1, Mesh::Point p2){
    return sqrt(square(p1[0]-p2[0]) + square(p1[1]-p2[1])+square(p1[2] - p2[2]));
}

Mesh::Scalar Geom_2D::square(Mesh::Scalar x){
    return x*x;
}

void Geom_2D::get_circumcircle(Point_2D &p0, Point_2D &p1, Point_2D &p2, double* pc, double &r){
    //void triangle_circumcircle ( double t[2*3], double &r, double pc[2] );

    double t[6];
    t[0] = (p0)[0];
    t[1] = (p0)[1];
    t[2] = (p1)[0];
    t[3] = (p1)[1];
    t[4] = (p2)[0];
    t[5] = (p2)[1];

    triangle_circumcircle(t, r, pc);
}

/**
 * @brief Geom_2D::get_2D_point
 * @param heh half edge in the mesh
 * @param p1 endpoint of the base line
 * @param p2 endpoint of the base line
 * @return the 2D equivalent to the mesh point, but flattened to the plane
 */
Point_2D Geom_2D::get_2D_point(Mesh::HalfedgeHandle heh, Point_2D &p1, Point_2D &p2){

    Mesh::Scalar angle = mesh->calc_sector_angle(heh);
    Mesh::Scalar len = mesh->calc_edge_length(heh);

    angle = - angle;
    //cout<<"angle: "<<angle<<endl;


    Mesh::HalfedgeLoopIter heIt;
    heIt = mesh->hl_begin(heh);
    heIt ++;

    Mesh::Scalar he_length = mesh->calc_edge_length(*heIt);

    //cout<<"next edge length: "<<he_length<<endl;


    double x = p1[0]-p2[0]; //make a vector of the edge we know, then rotate it.
    double y = p1[1]-p2[1];

    double xx = (x*cos(angle))-(y*sin(angle));  //apply rotation matrix
    double yy = (x*sin(angle))+(y*cos(angle));

    Point_2D p(xx,yy, 0.0);

    p.normalize();

    //cout<<"x: "<<x<<" y: "<<y<<" xx: "<<(*p)[0]<<" yy: "<<(*p)[1]<<endl;

    p *= he_length;

    p+=p2;

    //cout<<"x: "<<x<<" y: "<<y<<" xx: "<<(*p)[0]<<" yy: "<<(*p)[1]<<endl;

    return p;

}

vector<Point_2D> Geom_2D::test_flattening(Mesh::EdgeHandle ehandle){

    Mesh::HalfedgeHandle heh = mesh->halfedge_handle(ehandle, 0);

    //starting point. We position our edge along the x-axis.

    vector<Point_2D> p;

    cout<<"edge length: "<<mesh->calc_edge_length(heh)<<endl;

    p.push_back(Point_2D(0.0, 0.0, 0.0));  //0

    p.push_back(Point_2D(mesh->calc_edge_length(heh), 0.0, 0.0)); //1

    cout<<"point 2"<<endl;
    p.push_back(get_2D_point(heh, p[0], p[1])); //2

    Mesh::HalfedgeLoopIter heIt;
    heIt = mesh->hl_begin(heh);

    cout<<"point 3"<<endl;
    heIt ++;
    p.push_back(get_2D_point(mesh->opposite_halfedge_handle(*heIt), p[2], p[1]));  //3

    Mesh::HalfedgeHandle hedge= mesh->opposite_halfedge_handle(*heIt);
    Mesh::Point to, from;
    to = mesh->point(mesh->to_vertex_handle(hedge));
    from = mesh->point(mesh->from_vertex_handle(hedge));

    cout<<"base of point 3"<<endl;
    cout<<"point 1 x:"<<to[0]<<" y: "<<to[1]<<" z: "<<to[2]<<endl;
    cout<<"point 2 x:"<<from[0]<<" y: "<<from[1]<<" z: "<<from[2]<<endl;

    cout<<"point 4"<<endl;
    heIt ++;
    p.push_back(get_2D_point(mesh->opposite_halfedge_handle(*heIt), p[0], p[2])); //4

    Mesh::HalfedgeHandle heh2 = mesh->opposite_halfedge_handle(heh);

    cout<<"point 5"<<endl;

    p.push_back(get_2D_point(heh2, p[1], p[0])); //5
    heIt = mesh->hl_begin(heh2);

    cout<<"point 6"<<endl;
    heIt ++;
    p.push_back(get_2D_point(mesh->opposite_halfedge_handle(*heIt), p[5], p[0])); //6

    cout<<"point 7"<<endl;
    heIt ++;
    p.push_back(get_2D_point(mesh->opposite_halfedge_handle(*heIt), p[1], p[5])); //7

    for (int i = 0; i < 8; i++){
        cout<<(p[i])<<endl;
    }

    return p;
}



