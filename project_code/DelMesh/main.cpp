/*
 * delmesh.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: darryl
 */
#include "Delmesh.h"

void mesh_all(string input, string output);

int main(int argc, char* argv[]) {
    string input = "bunny.obj";
    string output = "output/ate-out.obj";

    //delmesh -i something -o something -ds stack -score 0 -nomedian
    bool out = false, all = false, median = true;
    int ds = 2, score = 2;
    int i = 1;
    while (argc > i){
        //input file
        cout<<"processing arguments: "<<argv[i]<<endl;
        if (strcmp(argv[i],"-i")==0){
            input = argv[++i];
            cout << "input: "<<input<<endl;
        }
        //output file
        else if (strcmp(argv[i],"-o")==0){
            out = true;
            output = argv[++i];
            cout << "output: "<< output<< endl;
        }
        //output all combinations of arguments for given input/output file
        else if (strcmp(argv[i],"-all")==0){
            all = true;
            cout<<"all"<<endl;
        }
        // choose data structure
        else if (strcmp(argv[i],"-ds")==0){
            i++;
            if(strcmp(argv[i],"queue")==0){
                ds = 0;}
            if(strcmp(argv[i],"pqueue")==0){
                ds = 1;
            }
            cout<<"ds: "<<ds<<endl;
        }
        //choose scoring system
        else if (strcmp(argv[i],"-s")==0){
            score = atoi(argv[++i]);
            cout<<"score: "<<score<<endl;

        }
        else if (strcmp(argv[i],"-nomedian")==0){
            cout<<"no median"<<endl;
            median = false;
        }else{
            cout<<"invalid argument"<<endl;
        }
        i++;
    }

    if (!out){
        output = input;
    }
    DelMesh* mesh;
    if (all){
        mesh_all(input, output);
        return 0;
    }

    mesh = new DelMesh(ds, score, input, output, median);
    mesh->process_mesh();
    delete(mesh);


}


//0 = priority queue
//1 = queue
//2 = stack
// 0 = score ++
// 1 = score +1.5
// 2 = score +5

//DelMesh mesh(DATA_STRUCTURE, score_type, input, "output/ate-out.obj");
//read the input mesh
//mesh.test_pq();
//mesh.test_face_colors();
//mesh.test_2D_flattening();
//mesh.test_samples();
//mesh.process_mesh();
//mesh.test_face_colors();

void mesh_all(string input, string output){

    int len = output.length();

    string out = output.substr(len-4,4);

    if (out.compare(".obj")==0){
        output = output.substr(0,len-4);
    }


    DelMesh* mesh = new DelMesh(0, 0, input, output + "-pq-0.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(1, 0, input, output + "-q-0.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(2, 0, input, output +"-s-0.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(0, 1, input, output +"-pq-1.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(1, 1, input, output +"-q-1.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(2, 1, input, output +"-s-1.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(0, 2, input, output +"-pq-2.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(1, 2, input, output +"-q-2.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(2, 2, input, output + "-s-2.obj", false);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(0, 0, input, output +"-pq-0-m.obj", true);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(1, 0, input, output +"-q-0-m.obj", true);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(2, 0, input, output +"-s-0-m.obj", true);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(0, 1, input, output +"-pq-1-m.obj", true);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(1, 1, input, output +"-q-1-m.obj", true);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(2, 1, input, output +"-s-1-m.obj", true);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(0, 2, input, output +"-pq-2-m.obj", true);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(1, 2, input, output +"-q-2-m.obj", true);
    mesh->process_mesh();
    delete(mesh);

    mesh = new DelMesh(2, 2, input, output + "-s-2-m.obj", true);
    mesh->process_mesh();
    delete(mesh);


}
