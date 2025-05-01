#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

#include "subdue.h"

void addEdge(Graph *graph, Parameters *parameters, ULONG source_id, ULONG target_id, BOOLEAN directed, std::string label_str)
{
    ULONG sourceVertexIndex;
    ULONG targetVertexIndex;

    // read and check vertex numbers
    if ((source_id < 1) || (source_id > graph->numVertices))
    {
        fprintf(stderr,
                "Error: reference to undefined vertex number %ld.\n",
                source_id);
        exit(1);
    }
    if ((target_id < 1) || (target_id > graph->numVertices))
    {
        fprintf(stderr,
                "Error: reference to undefined vertex number %ld.\n",
                target_id);
        exit(1);
    }
    sourceVertexIndex = source_id - 1;
    targetVertexIndex = target_id - 1;

    Label label;
    label.labelType = STRING_LABEL;
    label.labelValue.stringLabel = label_str.data();
    ULONG label_index = StoreLabel(&label, parameters->labelList);
    printf("Edge label: %s\n", parameters->labelList->labels[label_index].labelValue);

    Edge *newEdgeList;
    ULONG edgeListSize = graph->edgeListSize;

    // make sure there is enough room for another edge in the graph
    if (edgeListSize == graph->numEdges)
    {
        edgeListSize += LIST_SIZE_INC;
        newEdgeList = (Edge *)realloc(graph->edges, (sizeof(Edge) * (edgeListSize)));
        if (newEdgeList == NULL)
            OutOfMemoryError("AddEdge:newEdgeList");
        graph->edges = newEdgeList;
        graph->edgeListSize = edgeListSize;
    }

    // add edge to graph
    graph->edges[graph->numEdges].vertex1 = sourceVertexIndex;
    graph->edges[graph->numEdges].vertex2 = targetVertexIndex;
    graph->edges[graph->numEdges].label = label_index;
    graph->edges[graph->numEdges].directed = directed;
    graph->edges[graph->numEdges].used = FALSE;
    graph->edges[graph->numEdges].spansIncrement = FALSE;
    graph->edges[graph->numEdges].validPath = TRUE;

    AddEdgeToVertices(graph, graph->numEdges);

    ++graph->numEdges;
}

Graph *readSubdueGrpah(Parameters *parameters)
{
    Graph *graph;
    Vertex *newVertexList;
    ULONG *posEgsVertexIndices = NULL;
    ULONG numPosEgs = 0;

    numPosEgs++;
    posEgsVertexIndices = AddVertexIndex(posEgsVertexIndices, numPosEgs, 0);

    graph = AllocateGraph(0, 0);

    std::ifstream infile(parameters->inputFileName);
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        std::string entity;
        int id;
        ULONG source_id, target_id;

        iss >> entity;
        if (entity == "v")
        {
            double x = 0.0, y = 0.0, z = 0.0;
            std::string label_str;
            iss >> id >> label_str >> x >> y >> z;

            ULONG numVertices = graph->numVertices;
            ULONG vertexListSize = graph->vertexListSize;
            vertexListSize += LIST_SIZE_INC;
            newVertexList = (Vertex *)realloc(graph->vertices, (sizeof(Vertex) * (vertexListSize)));
            if (newVertexList == NULL)
                OutOfMemoryError("vertex list");
            graph->vertices = newVertexList;
            graph->vertexListSize = vertexListSize;

            Label label;
            label.labelType = STRING_LABEL;
            label.labelValue.stringLabel = label_str.data();
            label.level = 0;
            ULONG label_index = StoreLabel(&label, parameters->labelList);

            // store information in vertex
            graph->vertices[numVertices].label = label_index;
            graph->vertices[numVertices].numEdges = 0;
            graph->vertices[numVertices].edges = NULL;
            graph->vertices[numVertices].map = VERTEX_UNMAPPED;
            graph->vertices[numVertices].used = FALSE;
            graph->vertices[numVertices].position = gsl_vector_alloc(3);
            gsl_vector_set(graph->vertices[numVertices].position, 0, x);
            gsl_vector_set(graph->vertices[numVertices].position, 1, y);
            gsl_vector_set(graph->vertices[numVertices].position, 2, z);
            graph->vertices[numVertices].orientation = gsl_vector_alloc(3);
            gsl_vector_set(graph->vertices[numVertices].position, 0, 0);
            gsl_vector_set(graph->vertices[numVertices].position, 1, 0);
            gsl_vector_set(graph->vertices[numVertices].position, 2, 0);
            graph->vertices[numVertices].in_pattern = FALSE;
            ++graph->numVertices;
            std::cout << "Vertex: ID: " << id << " label: " << label_str << " position: ";
            printVector(graph->vertices[numVertices].position, 3);
            std::cout << "\n";
        }
        else if (entity == "u")
        {
            std::string label_str;
            iss >> source_id >> target_id >> label_str;
            addEdge(graph, parameters, source_id, target_id, FALSE, label_str);
            std::cout << "Undirected Edge: source: " << source_id << " target: " << target_id << " label: " << label_str << "\n";
        }
        else if (entity == "d")
        {
            std::string label_str;
            iss >> source_id >> target_id >> label_str;
            addEdge(graph, parameters, source_id, target_id, TRUE, label_str);
            std::cout << "Directed Edge: source: " << source_id << " target: " << target_id << " label: " << label_str << "\n";
        }
        else if (entity == "e")
        {
            std::string label_str;
            iss >> source_id >> target_id >> label_str;
            addEdge(graph, parameters, source_id, target_id, TRUE, label_str);
            std::cout << "Directed Edge: source: " << source_id << " target: " << target_id << " label: " << label_str << "\n";
        }
    }

    parameters->numPosEgs = numPosEgs;
    parameters->posEgsVertexIndices = posEgsVertexIndices;

    return graph;
}

int main(int argc, char *argv[])
{
    Parameters *parameters;
    parameters = GetParameters(argc, argv);
    parameters->posGraph = readSubdueGrpah(parameters);
    if (parameters->evalMethod == EVAL_MDL)
    {
        parameters->posGraphDL = MDL(parameters->posGraph,
                                     parameters->labelList->numLabels, parameters);
        if (parameters->negGraph != NULL)
        {
            parameters->negGraphDL =
                MDL(parameters->negGraph, parameters->labelList->numLabels,
                    parameters);
        }
    }

    PostProcessParameters(parameters);

    // Graph **compressed_graphs;
    // LabelList **all_label_lists;
    Graph *compressed_graphs[parameters->iterations];
    LabelList *all_label_lists[parameters->iterations];
    Substructure *discovered_subs[parameters->iterations];
    runSubdue(parameters, compressed_graphs, all_label_lists, discovered_subs);


    FreeParameters(parameters);

    return 0;
}
