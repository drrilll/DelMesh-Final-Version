/*
 * priorityqueue.h
 *
 *  Created on: Mar 11, 2016
 *      Author: darryl
 */

#ifndef PRIORITYQUEUE_H_
#define PRIORITYQUEUE_H_


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <iostream>
#include <queue>
#include "mytraits.h"


class custom_cmp{
public:
    custom_cmp(Mesh* mesh){
        this->mesh = mesh;
    }

    bool operator()(Mesh::EdgeHandle e1, Mesh::EdgeHandle e2){
        return (mesh->calc_edge_length(e1)<mesh->calc_edge_length(e2));
    }

private:
    Mesh* mesh;

};

/**
 * @brief The my_p_queue class
 * Depending on the integer argument given, this will act as a queue, a priority queue, or a stack
 */
class my_p_queue {
protected:
    std::deque<Mesh::EdgeHandle> c;
    custom_cmp* comp;

    /* 0: priority queue
     * 1: queue
     * 2: stack
     */
    int type = 0;

public:
    my_p_queue(Mesh* mesh, int type = 0)

{
        comp = new custom_cmp(mesh);
        if (type == 0){
            std::make_heap(c.begin(), c.end(), *comp);
        }
        this->type = type;
}

    ~my_p_queue(){
        delete comp;
    }

    void clear() {c.clear();}

    bool empty()       const { return c.empty();}

    std::size_t size() const { return c.size();}

    const Mesh::EdgeHandle& top()     const
    {
        if (type != 2){
            //queue or priority queue
            return c.front();
        }
        //it is a stack;
        return c.back();
    }

    void push(const Mesh::EdgeHandle& x)
    {
        c.push_back(x);
        if (type == 0){
            std::push_heap(c.begin(), c.end(), *comp);
        }
    }

    void pop()
    {
        //queue
        if (type == 1){
            c.pop_front();
            return;
        }
        if (type == 0){
            std::pop_heap(c.begin(), c.end(), *comp);
        }
        //stack or priority queue
        c.pop_back();
    }

};
#endif /* PRIORITYQUEUE_H_ */
