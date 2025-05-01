//---------------------------------------------------------------------------
// graphmatch.c
//
// Graph matching functions.
//
// Subdue 5
//---------------------------------------------------------------------------

#include "subdue.h"

#define HeapParent(i) ((i - 1) / 2)
#define HeapLeftChild(i) ((2 * i) + 1)
#define HeapRightChild(i) ((2 * i) + 2)

//---------------------------------------------------------------------------
// NAME:    GraphMatch
//
// INPUTS:  (Graph *g1)
//          (Graph *g2) - graphs to be matched
//          (LabelList *labelList) - list of vertex and edge labels
//          double threshold - upper bound on match cost
//          double *matchCost - pointer to pass back actual cost of match;
//                              ignored if NULL
//          (VertexMap *mapping) - array to hold final vertex mapping;
//                                 ignored if NULL
//
// RETURN:  (BOOLEAN) - TRUE is graphs match with cost less than threshold
//
// PURPOSE: Returns TRUE if g1 and g2 match with cost less than the given
// threshold.  If so, side-effects are to store the match cost in the
// variable pointed to by matchCost and to store the mapping between g1 and
// g2 in the given mapping input if non-NULL.
//---------------------------------------------------------------------------

BOOLEAN GraphMatch(Graph *g1, Graph *g2, LabelList *labelList,
                   double threshold, double *matchCost, VertexMap *mapping, BOOLEAN use_pose_cost)
{
    double cost;

    // first, quick check for exact matches
    if ((threshold == 0.0) &&
        ((g1->numVertices != g2->numVertices) ||
         (g1->numEdges != g2->numEdges)))
        return FALSE;

    // call InexactGraphMatch with larger graph first
    if (g1->numVertices < g2->numVertices)
        cost = InexactGraphMatch(g2, g1, labelList, threshold, mapping, use_pose_cost);
    else
        cost = InexactGraphMatch(g1, g2, labelList, threshold, mapping, use_pose_cost);

    // pass back actual match cost, if desired
    if (matchCost != NULL)
        *matchCost = cost;

    if (cost > threshold)
        return FALSE;

    return TRUE;
}

//---------------------------------------------------------------------------
// NAME:    InexactGraphMatch
//
// INPUTS:  Graph *g1
//          Graph *g2 - graphs to be matched
//          LabelList *labelList - list of vertex and edge labels
//          double threshold - upper bound on match cost
//          (VertexMap *mapping) - array to hold final vertex mapping;
//                                 if NULL, then ignored
//
// RETURN:  Cost of transforming g1 into an isomorphism of g2.  Will be
//          MAX_DOUBLE if cost exceeds threshold
//
// PURPOSE: Compute minimum-cost transformation of g1 to an
// isomorphism of g2, but any match cost exceeding the given threshold
// is not considerd.  Graph g1 should be the larger graph in terms of
// vertices.  A side-effect is to store the mapping between g1 and g2
// in the given mapping input if non-NULL.
//
// TODO: May want to input a partial mapping to influence mapping
// order of vertices in g1.
//---------------------------------------------------------------------------

double InexactGraphMatch(Graph *g1, Graph *g2, LabelList *labelList,
                         double threshold, VertexMap *mapping, BOOLEAN use_pose_cost)
{
    ULONG i, v1, v2;
    ULONG nv1 = g1->numVertices;
    ULONG nv2 = g2->numVertices;
    Edge *edge;
    MatchHeap *globalQueue = NULL;
    MatchHeap *localQueue = NULL;
    MatchHeapNode node;
    MatchHeapNode newNode;
    MatchHeapNode bestNode;
    double cost = 0.0;
    double newCost = 0.0;
    ULONG numNodes = 0;
    ULONG quickMatchThreshold = 0;
    BOOLEAN quickMatch = FALSE;
    BOOLEAN done = FALSE;
    ULONG *orderedVertices = NULL;
    ULONG *mapped1 = NULL; // mapping of vertices in g1 to vertices in g2
    ULONG *mapped2 = NULL; // mapping of vertices in g2 to vertices in g1

    // printf("Threshold: %f\n", threshold);

    // Compute threshold on mappings tried before changing from optimal
    // search to greedy search
    quickMatchThreshold = MaximumNodes(nv1);

    // Order vertices of g1 by degree
    orderedVertices = (ULONG *)malloc(sizeof(ULONG) * nv1);
    if (orderedVertices == NULL)
        OutOfMemoryError("orderedVertices");
    OrderVerticesByDegree(g1, orderedVertices);

    // Allocate arrays holding vertex maps from g1 to g2 and g2 to g1
    mapped1 = (ULONG *)malloc(sizeof(ULONG) * nv1);
    if (mapped1 == NULL)
        OutOfMemoryError("mapped1");
    mapped2 = (ULONG *)malloc(sizeof(ULONG) * nv2);
    if (mapped2 == NULL)
        OutOfMemoryError("mapped2");

    globalQueue = AllocateMatchHeap(nv1 * nv1);
    localQueue = AllocateMatchHeap(nv1);
    node.depth = 0;
    node.cost = 0.0;
    node.mapping = NULL;
    InsertMatchHeapNode(&node, globalQueue);
    bestNode.depth = 0;
    bestNode.cost = MAX_DOUBLE;
    bestNode.mapping = NULL;

    while ((!MatchHeapEmpty(globalQueue)) && (!done))
    {
        ExtractMatchHeapNode(globalQueue, &node);
        // if (node.depth >= nv1)
        // { // complete mapping found
        //     printf("node cost: %f, best cost: %f\n", node.cost, bestNode.cost);
        // }
        if (node.cost < bestNode.cost)
        {
            if (node.depth == nv1)
            { // complete mapping found
                // printf("node cost: %f, best cost: %f\n", node.cost, bestNode.cost);
                // printf("Complete mapping found\n");
                free(bestNode.mapping);
                bestNode.cost = node.cost;
                bestNode.depth = node.depth;
                bestNode.mapping = node.mapping;
                if (!quickMatch)
                    done = TRUE;
            }
            else
            { // expand node's partial mapping
              // compute mappings between graphs' vertices
                for (i = 0; i < nv1; i++)
                    mapped1[i] = VERTEX_UNMAPPED;
                for (i = 0; i < nv2; i++)
                    mapped2[i] = VERTEX_UNMAPPED;
                for (i = 0; i < node.depth; i++)
                {
                    mapped1[node.mapping[i].v1] = node.mapping[i].v2;
                    if (node.mapping[i].v2 != VERTEX_DELETED)
                        mapped2[node.mapping[i].v2] = node.mapping[i].v1;
                }
                v1 = orderedVertices[node.depth];
                // first, try mapping v1 to nothing
                double degree_multiplier = 1.0;
                if(g1->vertices[v1].degree >= 0)
                    degree_multiplier = 1.0 + g1->vertices[v1].degree / 6.0;  // The number in den = max degree in the graph
                // else
                //     degree_multiplier = 1.0 + label_wise_mean_degree[g1->vertices[v1].label]/6.0;
                // newCost = node.cost + DELETE_VERTEX_COST;
                newCost = node.cost + DELETE_VERTEX_COST * degree_multiplier;
                if ((newCost <= threshold) && (newCost < bestNode.cost))
                {
                    // printf("DELETE_VERTEX_COST survived\n");
                    for (i = 0; i < g1->vertices[v1].numEdges; i++)
                    {
                        edge = &g1->edges[g1->vertices[v1].edges[i]];
                        if (v1 == edge->vertex1)
                            v2 = edge->vertex2;
                        else
                            v2 = edge->vertex1;
                        if ((mapped1[v2] != VERTEX_DELETED) || (v2 == v1))
                        {
                            newCost += DELETE_EDGE_WITH_VERTEX_COST;
                        }
                    }
                }
                // if complete mapping, add cost for any unmapped vertices and
                // edges in g2
                if ((newCost <= threshold) && (newCost < bestNode.cost))
                {
                    // printf("DELETE_EDGE_WITH_VERTEX_COST survived\n");
                    if (node.depth == (nv1 - 1))
                    {
                        cost = InsertedVerticesCost(g2, mapped2);
                        newCost += cost;
                    }
                }
                if ((newCost <= threshold) && (newCost < bestNode.cost))
                {
                    // printf("InsertedVerticesCost survived\n");
                    // add new node to local queue
                    newNode.depth = node.depth + 1;
                    newNode.cost = newCost;
                    newNode.mapping =
                        AllocateNewMapping(node.depth + 1, node.mapping, v1,
                                           VERTEX_DELETED);
                    InsertMatchHeapNode(&newNode, localQueue);
                }
                // second, try mapping v1 to each unmapped vertex in g2
                for (v2 = 0; v2 < nv2; v2++)
                {
                    if (mapped2[v2] == VERTEX_UNMAPPED)
                    {
                        mapped1[v1] = v2;
                        mapped2[v2] = v1;
                        newCost = node.cost;
                        double pos_cost = 0.0;
                        double v1_mult = 1.0, v2_mult = 1.0, degree_mult = 1.0;
                        if(g1->vertices[v1].degree >= 0)
                            v1_mult = 1.0 + g1->vertices[v1].degree/6.0;
                        // else
                        //     v1_mult = label_wise_mean_degree[g1->vertices[v1].label];
                        
                        if(g2->vertices[v2].degree >= 0)
                            v2_mult = 1.0 + g2->vertices[v2].degree/6.0;
                        // else
                        //     v2_mult = label_wise_mean_degree[g2->vertices[v2].label];

                        if(v1_mult > v2_mult)
                            degree_multiplier = v1_mult;
                        else 
                            degree_mult = v2_mult;
                        newCost += SUBSTITUTE_VERTEX_LABEL_COST * degree_mult *
                                   LabelMatchFactor(g1->vertices[v1].label,
                                                    g2->vertices[v2].label, labelList);
                        
                        if ((newCost <= threshold) && (newCost < bestNode.cost))
                        {
                            cost = DeletedEdgesCost(g1, g2, v1, v2, mapped1, labelList);
                            newCost += cost;
                            cost = InsertedEdgesCost(g2, v2, mapped2);
                            newCost += cost;
                        }
                        // if complete mapping, add cost for any unmapped vertices
                        // and edges in g2
                        if ((newCost <= threshold) && (newCost < bestNode.cost))
                        {
                            if (node.depth == (nv1 - 1))
                            {
                                cost = InsertedVerticesCost(g2, mapped2);
                                newCost += cost;
                            }
                        }
                        // printf("Cost before pos cost: %g\n", newCost);
                        if (use_pose_cost && (newCost <= threshold) && (newCost < bestNode.cost))
                        {
                            // printf("v1: %lu, v2: %lu | node depth: %lu\n", v1, v2, node.depth);
                            pos_cost = PoseCost(g1, g2, v1, v2, mapped1, &node);
                            newCost += pos_cost;
                            // if(!strcmp(labelList->labels[g1->vertices[v1].label].labelValue.stringLabel, "SUB_1") || !strcmp(labelList->labels[g2->vertices[v2].label].labelValue.stringLabel, "SUB_1"))
                            // {
                            //     printf("Comparging\n");
                            //     PrintGraph(g1, labelList); PrintGraph(g2, labelList);
                            //     printf("Pos cost: %g \n", pos_cost);
                            // }
                        }

                        if ((newCost <= threshold) && (newCost < bestNode.cost))
                        {
                            // printf("UNMAPPED VERTICES COST survived\n");
                            // add new node to local queue
                            newNode.depth = node.depth + 1;
                            newNode.cost = newCost - pos_cost;
                            newNode.mapping =
                                AllocateNewMapping(node.depth + 1, node.mapping, v1, v2);
                            InsertMatchHeapNode(&newNode, localQueue);
                        }
                        mapped1[v1] = VERTEX_UNMAPPED;
                        mapped2[v2] = VERTEX_UNMAPPED;
                    }
                }
                free(node.mapping);
                // Add nodes in localQueue to globalQueue
                if (quickMatch)
                {
                    if (!MatchHeapEmpty(localQueue))
                    {
                        ExtractMatchHeapNode(localQueue, &node);
                        InsertMatchHeapNode(&node, globalQueue);
                        ClearMatchHeap(localQueue);
                    }
                }
                else
                    MergeMatchHeaps(localQueue, globalQueue); // clears localQueue
            }
        }
        else
            free(node.mapping);

        // check if maximum nodes exceeded, and if so, switch to greedy search
        numNodes++;
        if ((!quickMatch) && (numNodes > quickMatchThreshold))
        {
            CompressMatchHeap(globalQueue, nv1);
            quickMatch = TRUE;
        }
    } // end while

    // copy best mapping to input mapping array, if available
    if ((mapping != NULL) && (bestNode.mapping != NULL))
        for (i = 0; i < nv1; i++)
        {
            mapping[i].v1 = bestNode.mapping[i].v1;
            mapping[i].v2 = bestNode.mapping[i].v2;
        }
    double final_pos_cost = PoseCost(g1, g2, VERTEX_DELETED, VERTEX_DELETED, NULL, &bestNode);
    // printf("Pose cost: %f \n", final_pos_cost);
    if(use_pose_cost)
        bestNode.cost += final_pos_cost;
    // free memory
    free(bestNode.mapping);
    FreeMatchHeap(localQueue);
    FreeMatchHeap(globalQueue);
    free(mapped2);
    free(mapped1);
    free(orderedVertices);

    return bestNode.cost;
}

//---------------------------------------------------------------------------
// NAME: OrderVerticesByDegree
//
// INPUTS: (Graph *g) - graph whose vertices are to be sorted by degree
//         (ULONG *orderedVertices) - array to hold vertex indices
//                                    sorted by degree
//
// RETURN: (void)
//
// PURPOSE: Compute an array of vertex indices ordered by degree.
// Several studies have shown that matching higher degree vertices
// first will speed up the match.
//---------------------------------------------------------------------------

void OrderVerticesByDegree(Graph *g, ULONG *orderedVertices)
{
    ULONG nv = g->numVertices;
    ULONG i, j;
    ULONG degree;
    ULONG *vertexDegree = NULL;

    vertexDegree = (ULONG *)malloc(sizeof(ULONG) * nv);
    if (vertexDegree == NULL)
        OutOfMemoryError("vertexDegree");

    // insertion sort vertices by degree
    for (i = 0; i < nv; i++)
    {
        degree = g->vertices[i].numEdges;
        j = i;
        while ((j > 0) && (degree > vertexDegree[j - 1]))
        {
            vertexDegree[j] = vertexDegree[j - 1];
            orderedVertices[j] = orderedVertices[j - 1];
            j--;
        }
        vertexDegree[j] = degree;
        orderedVertices[j] = i;
    }
    free(vertexDegree);
}

//---------------------------------------------------------------------------
// NAME:    MaximumNodes
//
// INPUTS:  ULONG n
//
// RETURN:  n^MATCH_SEARCH_THRESHOLD_EXPONENT
//
// PURPOSE: Compute the maximum number of nodes (mappings) considered
// by InexactGraphMatch before switching from exhaustive to greedy
// search.  If MATCH_SEARCH_THRESHOLD_EXPONENT = 0.0, then assume
// exhaustive matching desired and retun maximum ULONG.
//---------------------------------------------------------------------------

ULONG MaximumNodes(ULONG n)
{
    if (MATCH_SEARCH_THRESHOLD_EXPONENT == 0.0)
        return MAX_UNSIGNED_LONG;
    else
        return (ULONG)pow((double)n, (double)MATCH_SEARCH_THRESHOLD_EXPONENT);
}

//---------------------------------------------------------------------------
// NAME: PoseCost
//
// INPUTS: (Graph *g1)
//         (Graph *g2) - graphs containing vertices being mapped
//         (ULONG v1) - vertex in g1 being mapped
//         (ULONG v2) - vertex in g2 being mapped
//         (ULONG *mapped1) - mapping of vertices in g1 to vertices in g2
//
// RETURN: (double) - cost of match edges according to given mapping
//
// PURPOSE: Compute the cost of matching the two vertices based on their
//          positions.
//---------------------------------------------------------------------------
double PoseCost(Graph *g1, Graph *g2, ULONG v1, ULONG v2,
                ULONG *mapped1, MatchHeapNode *node)
{
    const double match_distance_thr = 0.5;
    const double max_distance_error = 1.0;
    double cost = 0.0;

    VertexMap *filtered_map;
    ULONG filtered_depth = 0;
    /**/
    filtered_map = (VertexMap *)malloc(sizeof(VertexMap) * (node->depth+1));
    // printf("um: %lu, d: %lu\n", VERTEX_UNMAPPED, VERTEX_DELETED);
    // printf("Vertices being mapped: %lu, %lu\n", v1, v2);
    // printf("Node mappings:\n");
    for (int i = 0; i < node->depth; ++i)
    {
        if (node->mapping[i].v1 < VERTEX_DELETED && node->mapping[i].v2 < VERTEX_DELETED)
        {
            filtered_map[filtered_depth].v1 = node->mapping[i].v1;
            filtered_map[filtered_depth].v2 = node->mapping[i].v2;
            ++filtered_depth;
            // printf("%lu : %lu | ", node->mapping[i].v1, node->mapping[i].v2);
        }
    }

    if(v1 < VERTEX_DELETED && v2 < VERTEX_DELETED)
    {
        // filtered_map = (VertexMap *)realloc(filtered_map, filtered_depth * sizeof(VertexMap));
        filtered_map[filtered_depth].v1 = v1;
        filtered_map[filtered_depth].v2 = v2;
        ++filtered_depth;
        // printf("%lu : %lu | ", v1, v2);
    }
    /**/

    /**/
    // printf("\n");
    // // filtered_map = (VertexMap *)malloc(sizeof(VertexMap) * g1->numVertices);
    // for(ULONG v=0; v < g1->numVertices; ++v)
    // {
    //     if(mapped1[v] < VERTEX_UNMAPPED)
    //     {
    //         // filtered_map[filtered_depth].v1 = v;
    //         // filtered_map[filtered_depth].v2 = mapped1[v];
    //         // ++filtered_depth;
    //         printf("%lu : %lu | ", v, mapped1[v]);
    //     }
    // }
    // printf("\n");
    /**/


    // for (int i = 0; i < filtered_depth; ++i)
    // {
    //     printf("%d: Filtered map v1: %lu, v2: %lu\n", i, filtered_map[i].v1, filtered_map[i].v2);
    // }

    switch (filtered_depth)
    {
    case 0:
        cost = 0.0;
        break;
    case 1:
        cost = 0.0;
        break;
    case 2:
    {
        // printf("Case 1\n");
        // printf("v1: %ld, v2: %ld\n", node->mapping[0].v1, node->mapping[0].v2);
        gsl_vector *v0_1 = g1->vertices[filtered_map[0].v1].position;
        gsl_vector *v0_2 = g2->vertices[filtered_map[0].v2].position;
        gsl_vector *v1_1 = g1->vertices[filtered_map[1].v1].position;
        gsl_vector *v1_2 = g2->vertices[filtered_map[1].v2].position;
        gsl_vector *x1 = gsl_vector_alloc(3);
        gsl_vector *x2 = gsl_vector_alloc(3);
        // printf("0\n");
        // x1 = v1_1;
        // printf("1\n");
        gsl_vector_memcpy(x1, v1_1);
        // printf("0.5\n");
        // printf("x1[0]: %g\n", gsl_vector_get(x1,0));
        // printVector(x1, 1); printVector(v0_1, 3);
        gsl_vector_sub(x1, v0_1);
        // x2 = v1_2;
        // printf("2\n");
        gsl_vector_memcpy(x2, v1_2);
        gsl_vector_sub(x2, v0_2);
        // printf("1\n");
        double l1 = gsl_blas_dnrm2(x1);
        double l2 = gsl_blas_dnrm2(x2);
        // printf("2\n");
        double diff = abs(l1 - l2);
        if (diff <= match_distance_thr)
        {
            cost = 0.0;
        }
        else if (diff > match_distance_thr && diff <= max_distance_error)
        {
            cost = POSE_COST_MULTIPLIER * diff;
        }
        else
        {
            cost = max_distance_error * POSE_COST_MULTIPLIER;
        }
        // printf("Case 1: After calculations: g1: "); printVector(g1->vertices[v2].position, 3); printf(" g2: "); printVector(g2->vertices[v2].position, 3); printf("\n");
        // printf("3\n");
    }
    break;
    default:
    {
        // printf("Case default\n");
        gsl_matrix *m1 = gsl_matrix_alloc(filtered_depth, 3);
        gsl_matrix *m2 = gsl_matrix_alloc(filtered_depth, 3);
        gsl_matrix *m1_recentered = gsl_matrix_alloc(filtered_depth, 3);
        gsl_matrix *m2_recentered = gsl_matrix_alloc(filtered_depth, 3);

        for (ULONG i = 0; i < filtered_depth; ++i)
        {
            gsl_matrix_set_row(m1, i, g1->vertices[filtered_map[i].v1].position);
            gsl_matrix_set_row(m2, i, g2->vertices[filtered_map[i].v2].position);
            // printVector(g1->vertices[filtered_map[i].v1].position, 3); printf(" | "); printVector(g2->vertices[filtered_map[i].v2].position, 3); printf("\n");
        }

        // printf("Comparing subgraphs:\n");
        // printf("m1:\n");
        // printMatrix(m1, filtered_depth, 3);
        // printf("m2:\n");
        // printMatrix(m2, filtered_depth, 3);

        m1_recentered = recenterPoints(m1);
        m2_recentered = recenterPoints(m2);

        // printf("Recentered:\n");
        // printf("m1:\n");
        // printMatrix(m1_recentered, filtered_depth, 3);
        // printf("m2:\n");
        // printMatrix(m2_recentered, filtered_depth, 3);

        gsl_matrix *R = gsl_matrix_alloc(3, 3);
        R = getRotationMatrixSVD(m1_recentered, m2_recentered);
        // gsl_matrix_set_identity(R);
        // printf("R:\n");
        // printMatrix(R, R->size1, R->size2);

        gsl_matrix *m2_T = gsl_matrix_alloc(m2_recentered->size2, m2_recentered->size1);
        gsl_matrix_transpose_memcpy(m2_T, m2_recentered);
        gsl_matrix_transpose(R);
        gsl_matrix *m2_tfed = gsl_matrix_alloc(m2_recentered->size2, m2_recentered->size1);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, R, m2_T, 0.0, m2_tfed);
        gsl_matrix_transpose_memcpy(m2_recentered, m2_tfed);

        // printf("m2 tfed:\n");
        // printMatrix(m2_recentered, m2_recentered->size1, m2_recentered->size2);
        gsl_matrix_sub(m1_recentered, m2_recentered);

        // printf("Diff matrix:\n");
        // printMatrix(m1_recentered, m1_recentered->size1, m1_recentered->size2);

        double mean_diff = 0.0;
        gsl_vector *v = gsl_vector_alloc(3);
        for (int i = 0; i < m1_recentered->size1; ++i)
        {
            gsl_matrix_get_row(v, m1_recentered, i);
            mean_diff += gsl_blas_dnrm2(v);
        }
        mean_diff /= m1_recentered->size1;
        // printf("Mean diff %g\n", mean_diff);
        if (mean_diff <= match_distance_thr)
        {
            cost = 0.0;
        }
        else if (mean_diff > match_distance_thr && mean_diff <= max_distance_error)
        {
            cost = POSE_COST_MULTIPLIER * mean_diff;
        }
        else
        {
            cost = max_distance_error * POSE_COST_MULTIPLIER;
        }
        // cost = 0.0;

        gsl_matrix_free(m1_recentered);
        gsl_matrix_free(m2_recentered);
        // gsl_matrix_free(m1);
        // gsl_matrix_free(m2);
        // gsl_matrix_free(R);
        // gsl_matrix_free(m2_tfed);
        // gsl_matrix_free(m2_T);
    }
    break;
    }

    return cost;
}

// double PoseCost(Graph *g1, Graph *g2, ULONG v1, ULONG v2,
//                 ULONG *mapped1, MatchHeapNode *node)
// {
//     const double match_distance_thr = 0.5;
//     const double max_distance_error = 2.5;
//     double cost = 0.0;

//     VertexMap *filtered_map;
//     filtered_map = (VertexMap *)malloc(sizeof(VertexMap) * 3);
//     ULONG filtered_depth = 0;
//     for (int i = 0; i < node->depth; ++i)
//     {
//         if (filtered_depth >= 3)
//         {
//             break;
//         }
//         if (node->mapping[i].v1 >= 0)
//         {
//             if (node->mapping[i].v2 >= 0)
//             {
//                 ++filtered_depth;
//                 filtered_map[i].v1 = node->mapping[i].v1;
//                 filtered_map[i].v2 = node->mapping[i].v2;
//                 printf("%d: condition true v1: %d, v2: %d\n", i, node->mapping[i].v1, node->mapping[i].v2);
//             }
//         }
//     }

//     printf("filtered depth: %ld\n", filtered_depth);
//     for (int i = 0; i < filtered_depth; ++i)
//     {
//         printf("%d: v1: %ld, v2: %ld\n", i, node->mapping[i].v1, node->mapping[i].v2);
//     }

//     // printf("Did this print\n");
//     // printf("Node depth: %ld\n", node->depth);
//     switch (filtered_depth)
//     {
//     case 0:
//         cost = 0.0;
//         break;
//     case 1:
//     {
//         printf("Case 1\n");
//         // printf("v1: %ld, v2: %ld\n", node->mapping[0].v1, node->mapping[0].v2);
//         gsl_vector *v0_1 = g1->vertices[filtered_map[0].v1].position;
//         gsl_vector *v0_2 = g2->vertices[filtered_map[0].v2].position;
//         gsl_vector *v1_1 = g1->vertices[v1].position;
//         gsl_vector *v1_2 = g2->vertices[v2].position;
//         gsl_vector *x1 = gsl_vector_alloc(3);
//         gsl_vector *x2 = gsl_vector_alloc(3);
//         // printf("0\n");
//         // x1 = v1_1;
//         printf("1\n");
//         gsl_vector_memcpy(x1, v1_1);
//         // printf("0.5\n");
//         // printf("x1[0]: %g\n", gsl_vector_get(x1,0));
//         // printVector(x1, 1); printVector(v0_1, 3);
//         gsl_vector_sub(x1, v0_1);
//         // x2 = v1_2;
//         printf("2\n");
//         gsl_vector_memcpy(x2, v1_2);
//         gsl_vector_sub(x2, v0_2);
//         // printf("1\n");
//         double l1 = gsl_blas_dnrm2(x1);
//         double l2 = gsl_blas_dnrm2(x2);
//         // printf("2\n");
//         double diff = abs(l1 - l2);
//         if (diff <= match_distance_thr)
//         {
//             cost = 0.0;
//         }
//         else if (diff > match_distance_thr && diff <= max_distance_error)
//         {
//             cost = POSE_COST_MULTIPLIER * diff;
//         }
//         else
//         {
//             cost = max_distance_error * POSE_COST_MULTIPLIER;
//         }
//         // printf("Case 1: After calculations: g1: "); printVector(g1->vertices[v2].position, 3); printf(" g2: "); printVector(g2->vertices[v2].position, 3); printf("\n");
//         // printf("3\n");
//     }
//     break;
//     case 2:
//     {
//         printf("Case 2\n");
//         gsl_vector *v0_1 = g1->vertices[filtered_map[0].v1].position;
//         gsl_vector *v0_2 = g2->vertices[filtered_map[0].v2].position;
//         gsl_vector *v1_1 = g1->vertices[filtered_map[1].v1].position;
//         gsl_vector *v1_2 = g2->vertices[filtered_map[1].v2].position;
//         gsl_vector *v2_1 = g1->vertices[v1].position;
//         gsl_vector *v2_2 = g2->vertices[v2].position;

//         printf("mapping 0: v1: %d, v2: %d, mapping 1: v1: %d, v2: %d\n", filtered_map[0].v1, filtered_map[0].v2, filtered_map[1].v1, filtered_map[1].v2);
//         printf("Node label: %d\n", g2->vertices[filtered_map[1].v2].label);
//         printf("-1\n");
//         printf("Basis vectors 1:");
//         printVector(v0_1, 3);
//         printVector(v1_1, 3);
//         printVector(v2_1, 3);
//         printf("\n");
//         gsl_matrix *R1 = getRotationMatrixFromPlane(v0_1, v1_1, v2_1);
//         printf("0\n");
//         printVector(v0_2, 3);
//         // printVector(v1_2, 2);
//         printVector(v2_2, 3);
//         printf("\n");
//         printf("Basis vectors 2:");
//         printVector(v0_2, 3);
//         printVector(v1_2, 3);
//         printVector(v2_2, 3);
//         printf("\n");
//         gsl_matrix *R2 = getRotationMatrixFromPlane(v0_2, v1_2, v2_2);
//         printf("0.1======================\n");

//         gsl_vector *v2_1_tfed = gsl_vector_alloc(3);
//         gsl_vector *v2_2_tfed = gsl_vector_alloc(3);
//         gsl_blas_dgemv(CblasNoTrans, 1.0, R1, v2_1, 0.0, v2_1_tfed);
//         gsl_blas_dgemv(CblasNoTrans, 1.0, R2, v2_2, 0.0, v2_2_tfed);
//         printf("tfed vectors: \n");
//         printVector(v2_1_tfed, 3);
//         printVector(v2_2_tfed, 3);
//         printf("\n");

//         gsl_vector *x = gsl_vector_alloc(3);
//         // x = v2_1_tfed;
//         printf("1\n");
//         gsl_vector_memcpy(x, v2_1_tfed);
//         gsl_vector_sub(x, v2_2_tfed);
//         double diff = gsl_blas_dnrm2(x);
//         if (diff <= match_distance_thr)
//         {
//             cost = 0.0;
//         }
//         else if (diff > match_distance_thr && diff <= max_distance_error)
//         {
//             cost = POSE_COST_MULTIPLIER * diff;
//         }
//         else
//         {
//             cost = max_distance_error * POSE_COST_MULTIPLIER;
//         }

//         node->origin1 = gsl_vector_alloc(3);
//         node->R1 = gsl_matrix_alloc(3, 3);
//         printf("2\n");
//         gsl_vector_memcpy(node->origin1, v0_1);
//         copyMatrix(R1, node->R1, 3, 3);
//         // node->origin1 = v0_1;
//         // node->R1 = R1;
//         node->origin2 = gsl_vector_alloc(3);
//         node->R2 = gsl_matrix_alloc(3, 3);
//         printf("3\n");
//         gsl_vector_memcpy(node->origin2, v0_2);
//         copyMatrix(R2, node->R2, 3, 3);

//         // printf("Case 2: After calculations: g1: "); printVector(g1->vertices[v2].position, 3); printf(" g2: "); printVector(g2->vertices[v2].position, 3); printf("\n");
//         // node->origin2 = v0_2;
//         // node->R2 = R2;
//     }
//     break;
//     default:
//         if (node->depth <= 2)
//         {
//             printf("INCORRECT NODE DEPTH [%ld]", node->depth);
//         }
//         else
//         {
//             printf("Case 3\n");
//             gsl_vector *v2_1 = g1->vertices[v1].position;
//             printf("0.1\n");
//             gsl_vector *v2_2 = g2->vertices[v2].position;
//             printf("0.2\n");
//             gsl_vector *v2_1_tfed = gsl_vector_alloc(3);
//             printf("1\n");
//             gsl_vector_memcpy(v2_1_tfed, v2_1);
//             gsl_vector_sub(v2_1_tfed, node->origin1);
//             gsl_vector *v2_2_tfed = gsl_vector_alloc(3);
//             printf("2\n");
//             gsl_vector_memcpy(v2_2_tfed, v2_2);
//             gsl_vector_sub(v2_2_tfed, node->origin2);
//             gsl_blas_dgemv(CblasNoTrans, 1.0, node->R1, v2_1, 0.0, v2_1_tfed);
//             gsl_blas_dgemv(CblasNoTrans, 1.0, node->R2, v2_2, 0.0, v2_2_tfed);

//             gsl_vector *x = gsl_vector_alloc(3);
//             // x = v2_1_tfed;
//             printf("3\n");
//             gsl_vector_memcpy(x, v2_1_tfed);
//             gsl_vector_sub(x, v2_2_tfed);
//             double diff = gsl_blas_dnrm2(x);
//             if (diff <= match_distance_thr)
//             {
//                 cost = 0.0;
//             }
//             else if (diff > match_distance_thr && diff <= max_distance_error)
//             {
//                 cost = POSE_COST_MULTIPLIER * diff;
//             }
//             else
//             {
//                 cost = max_distance_error * POSE_COST_MULTIPLIER;
//             }

//             printf("Case 3: After calculations: g1: ");
//             printVector(g1->vertices[v2].position, 3);
//             printf(" g2: ");
//             printVector(g2->vertices[v2].position, 3);
//             printf("\n");
//         }
//         break;
//     }

//     return cost;
// }

//---------------------------------------------------------------------------
// NAME: DeletedEdgesCost
//
// INPUTS: (Graph *g1)
//         (Graph *g2) - graphs containing vertices being mapped
//         (ULONG v1) - vertex in g1 being mapped
//         (ULONG v2) - vertex in g2 being mapped
//         (ULONG *mapped1) - mapping of vertices in g1 to vertices in g2
//         (LabelList *labelList) - label list containing labels for g1 and g2
//
// RETURN: (double) - cost of match edges according to given mapping
//
// PURPOSE: Compute the cost of matching edges involved in the new
// mapping, which has just added v1 -> v2.  In the case of multiple
// edges between two vertices, do a greedy search to find a low-cost
// mapping of edges to edges.
//
// NOTE: Assumes InsertedEdgesCost() run right after this one.
//---------------------------------------------------------------------------

double DeletedEdgesCost(Graph *g1, Graph *g2, ULONG v1, ULONG v2,
                        ULONG *mapped1, LabelList *labelList)
{
    ULONG e1, e2;
    Edge *edge1, *edge2;
    ULONG otherVertex1, otherVertex2;
    Edge *bestMatchEdge;
    double bestMatchCost;
    double matchCost;
    double totalCost = 0.0;

    // try to match each edge involving v1 to an edge involving v2
    for (e1 = 0; e1 < g1->vertices[v1].numEdges; e1++)
    {
        edge1 = &g1->edges[g1->vertices[v1].edges[e1]];
        if (edge1->vertex1 == v1)
            otherVertex1 = edge1->vertex2;
        else
            otherVertex1 = edge1->vertex1;
        if ((mapped1[otherVertex1] != VERTEX_UNMAPPED) &&
            (mapped1[otherVertex1] != VERTEX_DELETED))
        {
            // target vertex of edge also mapped
            bestMatchEdge = NULL;
            bestMatchCost = -1.0;
            otherVertex2 = mapped1[otherVertex1];
            for (e2 = 0; e2 < g2->vertices[v2].numEdges; e2++)
            {
                edge2 = &g2->edges[g2->vertices[v2].edges[e2]];
                if ((!edge2->used) &&
                    (((edge2->vertex1 == otherVertex2) && (edge2->vertex2 == v2)) ||
                     ((edge2->vertex1 == v2) && (edge2->vertex2 == otherVertex2))))
                {
                    // compute cost of matching edge in g1 to edge in g2
                    matchCost = 0.0;
                    if (edge1->directed != edge2->directed)
                        matchCost += SUBSTITUTE_EDGE_DIRECTION_COST;
                    if (edge1->directed && edge2->directed &&
                        (edge1->vertex1 != edge1->vertex2) &&
                        (((edge1->vertex1 == v1) &&
                          (edge2->vertex1 == otherVertex2)) ||
                         ((edge1->vertex1 == otherVertex1) &&
                          (edge2->vertex1 == v2))))
                        matchCost += REVERSE_EDGE_DIRECTION_COST;
                    matchCost += SUBSTITUTE_EDGE_LABEL_COST *
                                 LabelMatchFactor(edge1->label, edge2->label, labelList);
                    if ((matchCost < bestMatchCost) || (bestMatchCost < 0.0))
                    {
                        bestMatchCost = matchCost;
                        bestMatchEdge = edge2;
                    }
                }
            }
            // if matching edge found, then add cost of match and mark edge used;
            // else add cost of deleting edge from g1
            if (bestMatchEdge != NULL)
            {
                bestMatchEdge->used = TRUE;
                totalCost += bestMatchCost;
            }
            else
                totalCost += DELETE_EDGE_COST;
        }
    }
    return totalCost;
}

//---------------------------------------------------------------------------
// NAME: InsertedEdgesCost
//
// INPUTS: (Graph *g2) - graph containing vertex being mapped to
//         (ULONG v2) - vertex in g2 being mapped to
//         (ULONG *mapped2) - array mapping vertices of g2 to vertices of g1
//
// RETURN: (double) - cost of inserting edges found in g2 between v2
//                    and another mapped vertex, but not matched to
//                    edges in g1
//
// PURPOSE: Compute the cost of adding edges to the transformation
// from g1 to g2 for edges in g2 that are between v2 and another mapped
// vertex, but are not matched to edges in g1 by DeleteEdgesCost().
//
// NOTE: Assumes DeletedEdgesCost() run before this one.
//---------------------------------------------------------------------------

double InsertedEdgesCost(Graph *g2, ULONG v2, ULONG *mapped2)
{
    ULONG e2;
    Edge *edge2;
    double totalCost = 0.0;

    for (e2 = 0; e2 < g2->vertices[v2].numEdges; e2++)
    {
        edge2 = &g2->edges[g2->vertices[v2].edges[e2]];
        if ((!edge2->used) &&
            (mapped2[edge2->vertex1] != VERTEX_UNMAPPED) &&
            (mapped2[edge2->vertex2] != VERTEX_UNMAPPED))
            totalCost += INSERT_EDGE_COST;
        edge2->used = FALSE;
    }
    return totalCost;
}

//---------------------------------------------------------------------------
// NAME: InsertedVerticesCost
//
// INPUTS: (Graph *g2) - graph being transformed to that may contain vertices
//                       not mapped to
//         (ULONG *mapped2) - mapping of vertices from g2 to g1
//
// RETURN: (double) - cost of inserting vertices and edges not mapped in g2
//
// PURPOSE: When all vertices in g1 have been mapped, some of them may have
// been mapped to NODE_DELETED and some of g2's vertices may then be
// unmapped.  This function computes the cost of inserting the necessary
// vertices and edges to transform g1 into g2.  Care is taken not to
// inflate the cost for edges both of whose vertices are unmapped.
//---------------------------------------------------------------------------

double InsertedVerticesCost(Graph *g2, ULONG *mapped2)
{
    ULONG v;
    ULONG v2;
    ULONG e;
    double cost = 0.0;

    for (v = 0; v < g2->numVertices; v++)
        if (mapped2[v] == VERTEX_UNMAPPED)
        {
            cost += INSERT_VERTEX_COST;
            // add cost for edges to mapped vertices and self edges
            for (e = 0; e < g2->vertices[v].numEdges; e++)
            {
                if (v == g2->edges[g2->vertices[v].edges[e]].vertex1)
                    v2 = g2->edges[g2->vertices[v].edges[e]].vertex2;
                else
                    v2 = g2->edges[g2->vertices[v].edges[e]].vertex1;
                if ((mapped2[v2] != VERTEX_UNMAPPED) || (v2 == v))
                    cost += INSERT_EDGE_WITH_VERTEX_COST;
            }
        }
    return cost;
}

//---------------------------------------------------------------------------
// Match Node Heap Functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// NAME: AllocateMatchHeap
//
// INPUTS: (ULONG size) - nodes of memory initially allocated to heap
//
// RETURN: (MatchHeap *) - pointer to newly-allocated empty match node heap
//
// PURPOSE: Allocate memory for a match node heap with enough room for size
// nodes.
//---------------------------------------------------------------------------

MatchHeap *AllocateMatchHeap(ULONG size)
{
    MatchHeap *heap;

    heap = (MatchHeap *)malloc(sizeof(MatchHeap));
    if (heap == NULL)
        OutOfMemoryError("AllocateMatchHeap:MatchHeap");
    heap->numNodes = 0;
    heap->size = size;
    heap->nodes = (MatchHeapNode *)malloc(size * sizeof(MatchHeapNode));
    if (heap->nodes == NULL)
        OutOfMemoryError("AllocateMatchHeap:heap->nodes");

    return heap;
}

//---------------------------------------------------------------------------
// NAME: AllocateNewMapping
//
// INPUTS: (ULONG depth) - depth (vertices mapped) of new node in search space
//         (VertexMap *mapping) - partial mapping to be augmented
//         (ULONG v1)
//         (ULONG v2) - node's mapping includes v1 -> v2
//
// RETURN:  (VertexMap *) - newly-allocated vertex mapping
//
// PURPOSE: Allocate a new vertex mapping as a copy of the given mapping
// appended with the new map v1 -> v2.
//---------------------------------------------------------------------------

VertexMap *AllocateNewMapping(ULONG depth, VertexMap *mapping,
                              ULONG v1, ULONG v2)
{
    VertexMap *newMapping;
    ULONG i;

    newMapping = (VertexMap *)malloc(sizeof(VertexMap) * depth);
    if (newMapping == NULL)
        OutOfMemoryError("AllocateNewMapping: newMapping");
    for (i = 0; i < (depth - 1); i++)
    {
        newMapping[i].v1 = mapping[i].v1;
        newMapping[i].v2 = mapping[i].v2;
    }
    newMapping[depth - 1].v1 = v1;
    newMapping[depth - 1].v2 = v2;

    return newMapping;
}

//---------------------------------------------------------------------------
// NAME: MatchHeapEmpty
//
// INPUTS: (MatchHeap *heap) - heap to check if empty
//
// RETURN: (BOOLEAN) - TRUE if heap is empty; else FALSE
//
// PURPOSE: Determine if given match node heap is empty.
//---------------------------------------------------------------------------

BOOLEAN MatchHeapEmpty(MatchHeap *heap)
{
    if (heap->numNodes == 0)
        return TRUE;
    return FALSE;
}

//---------------------------------------------------------------------------
// NAME: InsertMatchHeapNode
//
// INPUTS: (MatchHeapNode *node) - node to insert
//         (MatchHeap *heap) - heap to insert into
//
// RETURN:  (void)
//
// PURPOSE: Insert given node into given heap, maintaining increasing
// order by cost, and for nodes with the same cost, in decreasing order
// by depth.
//---------------------------------------------------------------------------

void InsertMatchHeapNode(MatchHeapNode *node, MatchHeap *heap)
{
    ULONG i;
    ULONG parent;
    MatchHeapNode *node2;
    BOOLEAN done;

    heap->numNodes++;
    // add more memory to heap if necessary
    if (heap->numNodes > heap->size)
    {
        heap->size = 2 * heap->size;
        heap->nodes = (MatchHeapNode *)realloc(heap->nodes, heap->size * sizeof(MatchHeapNode));
        if (heap->nodes == NULL)
            OutOfMemoryError("InsertMatchHeapNode:heap->nodes");
    }
    // find right place for new node
    i = heap->numNodes - 1;
    done = FALSE;
    while ((i > 0) && (!done))
    {
        parent = HeapParent(i);
        node2 = &heap->nodes[parent];
        if ((node->cost < node2->cost) ||
            ((node->cost == node2->cost) && (node->depth > node2->depth)))
        {
            heap->nodes[i].cost = heap->nodes[parent].cost;
            heap->nodes[i].depth = heap->nodes[parent].depth;
            heap->nodes[i].mapping = heap->nodes[parent].mapping;
            i = parent;
        }
        else
            done = TRUE;
    }
    // store new node
    heap->nodes[i].cost = node->cost;
    heap->nodes[i].depth = node->depth;
    heap->nodes[i].mapping = node->mapping;
}

//---------------------------------------------------------------------------
// NAME: ExtractMatchHeapNode
//
// INPUTS: (MatchHeap *heap) - heap to be removed from
//         (MatchHeapNode *node) - node into which to copy information
//                                 from best node
//
// RETURN: (void)
//
// PURPOSE: Return, in the given node, the best node in the heap and
// maintain ordering.
//---------------------------------------------------------------------------

void ExtractMatchHeapNode(MatchHeap *heap, MatchHeapNode *node)
{
    ULONG i;

    // copy best node to input storage node
    node->cost = heap->nodes[0].cost;
    node->depth = heap->nodes[0].depth;
    node->mapping = heap->nodes[0].mapping;

    // copy last node in heap array to first
    i = heap->numNodes - 1;
    heap->nodes[0].cost = heap->nodes[i].cost;
    heap->nodes[0].depth = heap->nodes[i].depth;
    heap->nodes[0].mapping = heap->nodes[i].mapping;
    heap->numNodes--;

    HeapifyMatchHeap(heap);
}

//---------------------------------------------------------------------------
// NAME: HeapifyMatchHeap
//
// INPUTS: (MatchHeap *heap) - heap to heapify from root
//
// RETURN:  (void)
//
// PURPOSE: Restores the heap property of the heap starting at the root
// node.  The heap property is that parent nodes have less cost than their
// children, and if the cost is the same, then parents have greater or
// equal depth than their children.
//---------------------------------------------------------------------------

void HeapifyMatchHeap(MatchHeap *heap)
{
    ULONG parent;
    ULONG best;
    ULONG leftChild;
    ULONG rightChild;
    MatchHeapNode *parentNode;
    MatchHeapNode *bestNode;
    MatchHeapNode *leftNode;
    MatchHeapNode *rightNode;
    ULONG tmpDepth;
    double tmpCost;
    VertexMap *tmpMapping;

    parent = 0;
    best = 1;
    while (parent != best)
    {
        parentNode = &heap->nodes[parent];
        leftChild = HeapLeftChild(parent);
        rightChild = HeapRightChild(parent);
        best = parent;
        bestNode = parentNode;
        // check if left child better than parent
        if (leftChild < heap->numNodes)
        {
            leftNode = &heap->nodes[leftChild];
            if ((leftNode->cost < parentNode->cost) ||
                ((leftNode->cost == parentNode->cost) &&
                 (leftNode->depth > parentNode->depth)))
            {
                best = leftChild;
                bestNode = leftNode;
            }
        }
        // check if right child better than best so far
        if (rightChild < heap->numNodes)
        {
            rightNode = &heap->nodes[rightChild];
            if ((rightNode->cost < bestNode->cost) ||
                ((rightNode->cost == bestNode->cost) &&
                 (rightNode->depth > bestNode->depth)))
            {
                best = rightChild;
                bestNode = rightNode;
            }
        }
        // if child better than parent, then swap
        if (parent != best)
        {
            tmpCost = parentNode->cost;
            tmpDepth = parentNode->depth;
            tmpMapping = parentNode->mapping;
            parentNode->cost = bestNode->cost;
            parentNode->depth = bestNode->depth;
            parentNode->mapping = bestNode->mapping;
            bestNode->cost = tmpCost;
            bestNode->depth = tmpDepth;
            bestNode->mapping = tmpMapping;
            parent = best;
            best = 0; // something other than parent so while loop continues
        }
    }
}

//---------------------------------------------------------------------------
// NAME: MergeMatchHeaps
//
// INPUTS: (MatchHeap *heap1) - match node list to be inserted from
//         (MatchHeap *heap2) - match node list to be inserted into
//
// RETURN:  (void)
//
// PURPOSE: Insert the match nodes of heap1 into heap2 while maintaining
// the order of heap2.  Ordering is by increasing cost, or if costs are
// equal, by decreasing depth.  The nodes of heap2 are removed, but the
// memory remains allocated.
//---------------------------------------------------------------------------

void MergeMatchHeaps(MatchHeap *heap1, MatchHeap *heap2)
{
    ULONG i;
    MatchHeapNode *node;

    for (i = 0; i < heap1->numNodes; i++)
    {
        node = &heap1->nodes[i];
        InsertMatchHeapNode(node, heap2);
        node->mapping = NULL;
    }
    heap1->numNodes = 0;
}

//---------------------------------------------------------------------------
// NAME: CompressMatchHeap
//
// INPUTS: (MatchHeap *heap) - match node heap to compress
//         (ULONG n) - limit used for bounding match node heap
//
// RETURN: (void)
//
// PURPOSE: Compress match node heap for the beginning of greedy
// search within the InexactGraphMatch function.  The first n nodes
// are left on the heap.  If there are more nodes on the heap, then
// the nodes with unique costs remain on the heap, and the rest are
// freed.  Note that the heap is assumed to already be in increasing
// order by cost, and for nodes having the same cost, in decreasing
// order by depth.
//---------------------------------------------------------------------------

void CompressMatchHeap(MatchHeap *heap, ULONG n)
{
    MatchHeap *newHeap;
    MatchHeapNode node1;
    MatchHeapNode node2;

    newHeap = AllocateMatchHeap(n);

    // keep best n nodes
    while ((n > 0) && (!MatchHeapEmpty(heap)))
    {
        ExtractMatchHeapNode(heap, &node1);
        InsertMatchHeapNode(&node1, newHeap);
        n--;
    }

    // keep remaining nodes with unique costs
    while (!MatchHeapEmpty(heap))
    {
        ExtractMatchHeapNode(heap, &node2);
        if (node1.cost == node2.cost)
            free(node2.mapping);
        else
        {
            InsertMatchHeapNode(&node2, newHeap);
            node1.cost = node2.cost;
        }
    }

    // free old heap and install new one in its place
    free(heap->nodes);
    heap->nodes = newHeap->nodes;
    heap->size = newHeap->size;
    heap->numNodes = newHeap->numNodes;
    free(newHeap);
}

//---------------------------------------------------------------------------
// NAME: PrintMatchHeapNode
//
// INPUTS: (MatchHeapNode *node) - match node to print
//
// RETURN: (void)
//
// PURPOSE: Print match node.
//---------------------------------------------------------------------------

void PrintMatchHeapNode(MatchHeapNode *node)
{
    ULONG i;
    printf("MatchHeapNode: depth = %lu, cost = %f, mapping =",
           node->depth, node->cost);
    if (node->depth > 0)
    {
        printf("\n");
        for (i = 0; i < node->depth; i++)
        {
            printf("            %lu -> ", node->mapping[i].v1);
            if (node->mapping[i].v2 == VERTEX_UNMAPPED)
                printf("unmapped\n");
            else if (node->mapping[i].v2 == VERTEX_DELETED)
                printf("deleted\n");
            else
                printf("%lu\n", node->mapping[i].v2);
        }
    }
    else
        printf(" NULL\n");
}

//---------------------------------------------------------------------------
// NAME: PrintMatchHeap
//
// INPUTS: (MatchHeap *heap) - match node heap to print
//
// RETURN: (void)
//
// PURPOSE: Print match node list.
//---------------------------------------------------------------------------

void PrintMatchHeap(MatchHeap *heap)
{
    ULONG i;
    MatchHeapNode *node;

    printf("MatchHeap:\n");
    for (i = 0; i < heap->numNodes; i++)
    {
        node = &heap->nodes[i];
        printf("(%lu) ", i);
        PrintMatchHeapNode(node);
    }
}

//---------------------------------------------------------------------------
// NAME: ClearMatchHeap
//
// INPUTS: (MatchHeap *heap) - match node heap to clear
//
// RETURN: (void)
//
// PURPOSE: Free mappings of match nodes still in heap and reset heap to
// have zero nodes.
//---------------------------------------------------------------------------

void ClearMatchHeap(MatchHeap *heap)
{
    ULONG i;

    for (i = 0; i < heap->numNodes; i++)
        free(heap->nodes[i].mapping);
    heap->numNodes = 0;
}

//---------------------------------------------------------------------------
// NAME: FreeMatchHeap
//
// INPUTS: (MatchHeap *heap) - match node heap to free
//
// RETURN: (void)
//
// PURPOSE: Free memory in given match node heap.
//---------------------------------------------------------------------------

void FreeMatchHeap(MatchHeap *heap)
{
    ClearMatchHeap(heap);
    free(heap->nodes);
    free(heap);
}

void crossProduct(gsl_vector *v1, gsl_vector *v2, gsl_vector *v_out)
{
    gsl_vector_set(v_out, 0, (gsl_vector_get(v1, 1) * gsl_vector_get(v2, 2) - gsl_vector_get(v1, 2) * gsl_vector_get(v2, 1)));
    gsl_vector_set(v_out, 1, (gsl_vector_get(v1, 2) * gsl_vector_get(v2, 0) - gsl_vector_get(v1, 0) * gsl_vector_get(v2, 2)));
    gsl_vector_set(v_out, 2, (gsl_vector_get(v1, 0) * gsl_vector_get(v2, 1) - gsl_vector_get(v1, 1) * gsl_vector_get(v2, 0)));
}

void printVector(gsl_vector *v, int num_elems)
{
    printf("[");
    for (int i = 0; i < num_elems; i++)
    {
        printf("%g ", gsl_vector_get(v, i));
    }
    printf("]");
}

void printMatrix(gsl_matrix *m, int num_rows, int num_cols)
{
    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            printf("%g ", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }
}

void normalize(gsl_vector *v)
{
    gsl_vector_scale(v, (1.0 / gsl_blas_dnrm2(v)));
}

void copyMatrix(gsl_matrix *src, gsl_matrix *dst, int num_ros, int num_cols)
{
    gsl_vector *r = gsl_vector_alloc(num_cols);
    for (int i = 0; i < num_ros; ++i)
    {
        gsl_matrix_get_row(r, src, i);
        gsl_matrix_set_row(dst, i, r);
    }
}

void transformToPlaneOf(gsl_vector *v0, gsl_vector *v1, gsl_vector *v2, gsl_vector *v_in, gsl_vector *v_out)
{
    gsl_vector *x = gsl_vector_alloc(3);
    gsl_vector *y = gsl_vector_alloc(3);
    gsl_vector *z = gsl_vector_alloc(3);
    gsl_vector *vx = gsl_vector_alloc(3);
    gsl_vector *vy = gsl_vector_alloc(3);
    gsl_vector *vz = gsl_vector_alloc(3);
    gsl_vector_set_basis(vx, 0);
    gsl_vector_set_basis(vy, 1);
    gsl_vector_set_basis(vz, 2);

    // x = v1 - v0
    x = v1;
    gsl_vector_sub(x, v0);

    gsl_vector *v_1_0 = v1;
    gsl_vector_sub(v_1_0, v0);
    gsl_vector *v_2_0 = v2;
    gsl_vector_sub(v_2_0, v0);

    // z = v_1_0 (x) v_2_0
    crossProduct(v_1_0, v_2_0, z);
    // y = z (x) x
    crossProduct(z, x, y);

    normalize(x);
    normalize(y);
    normalize(z);
    // printf("x: ");
    // printVector(x, 3);
    // printf(" y: ");
    // printVector(y, 3);
    // printf(" z: ");
    // printVector(z, 3);
    // printf("\n");
    // printf("vx: ");
    // printVector(vx, 3);
    // printf(" vy: ");
    // printVector(vy, 3);
    // printf(" vz: ");
    // printVector(vz, 3);
    // printf(" \n");

    double dot_prod;
    gsl_blas_ddot(x, v_in, &dot_prod);
    gsl_vector_set(v_out, 0, dot_prod);
    gsl_blas_ddot(y, v_in, &dot_prod);
    gsl_vector_set(v_out, 1, dot_prod);
    gsl_blas_ddot(z, v_in, &dot_prod);
    gsl_vector_set(v_out, 2, dot_prod);
}

gsl_matrix *getRotationMatrixFromPlane(gsl_vector *v0, gsl_vector *v1, gsl_vector *v2)
{
    // gsl_vector *a = gsl_vector_alloc(3);
    gsl_vector *x = gsl_vector_alloc(3);
    gsl_vector *y = gsl_vector_alloc(3);
    gsl_vector *z = gsl_vector_alloc(3);
    gsl_vector *vx = gsl_vector_alloc(3);
    gsl_vector *vy = gsl_vector_alloc(3);
    gsl_vector *vz = gsl_vector_alloc(3);
    gsl_vector_set_basis(vx, 0);
    gsl_vector_set_basis(vy, 1);
    gsl_vector_set_basis(vz, 2);

    // x = v1 - v0
    // x = v1;
    x = gsl_vector_alloc(3);
    gsl_vector_memcpy(x, v1);
    gsl_vector_sub(x, v0);

    // gsl_vector *v_1_0 = gsl_vector_alloc(3);
    // gsl_vector_memcpy(v_1_0, v1);
    // gsl_vector_sub(v_1_0, v0);
    gsl_vector *v_2_0 = gsl_vector_alloc(3);
    gsl_vector_memcpy(v_2_0, v2);
    gsl_vector_sub(v_2_0, v0);

    // z = v_1_0 (x) v_2_0
    crossProduct(x, v_2_0, z);
    // y = z (x) x
    crossProduct(z, x, y);

    // printf("Before normalization: x: ");
    // printVector(x, 3);
    // printf(" y: ");
    // printVector(y, 3);
    // printf(" z: ");
    // printVector(z, 3);
    // printf("\n");
    normalize(x);
    normalize(y);
    normalize(z);
    // printf("x: ");
    // printVector(x, 3);
    // printf(" y: ");
    // printVector(y, 3);
    // printf(" z: ");
    // printVector(z, 3);
    // printf("\n");
    // printf("vx: ");
    // printVector(vx, 3);
    // printf(" vy: ");
    // printVector(vy, 3);
    // printf(" vz: ");
    // printVector(vz, 3);
    // printf(" \n");

    // double dot_prod;
    // gsl_blas_ddot(x, v_in, &dot_prod);
    // gsl_vector_set(v_out, 0, dot_prod);
    // gsl_blas_ddot(y, v_in, &dot_prod);
    // gsl_vector_set(v_out, 1, dot_prod);
    // gsl_blas_ddot(z, v_in, &dot_prod);
    // gsl_vector_set(v_out, 2, dot_prod);

    gsl_matrix *R = gsl_matrix_alloc(3, 3);
    gsl_matrix_set_row(R, 0, x);
    gsl_matrix_set_row(R, 1, y);
    gsl_matrix_set_row(R, 2, z);

    return R;
}

gsl_matrix *recenterPoints(gsl_matrix *m)
{
    gsl_vector *vr = gsl_vector_alloc(m->size2);
    gsl_vector *vr2 = gsl_vector_alloc(m->size2);
    gsl_matrix_get_row(vr, m, 0);
    for (int i = 1; i < m->size1; ++i)
    {
        gsl_matrix_get_row(vr2, m, i);
        gsl_vector_add(vr, vr2);
    }

    gsl_vector_scale(vr, 1.0 / (double)m->size1); // Calculating centroid
    gsl_matrix *msub = gsl_matrix_alloc(m->size1, m->size2);
    for (int i = 0; i < m->size1; ++i)
    {
        gsl_matrix_set_row(msub, i, vr);
    }

    gsl_matrix *m_recentered = gsl_matrix_alloc(m->size1, m->size2);
    gsl_matrix_memcpy(m_recentered, m);
    gsl_matrix_sub(m_recentered, msub);

    return m_recentered;
}

gsl_matrix *getRotationMatrixSVD(gsl_matrix *m1, gsl_matrix *m2)
{
    gsl_matrix *m1_T = gsl_matrix_alloc(m1->size2, m1->size1);
    // gsl_matrix_transpose(m1);
    gsl_matrix_transpose_memcpy(m1_T, m1);
    gsl_matrix *H = gsl_matrix_alloc(3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m1_T, m2, 0.0, H);

    gsl_matrix *U = gsl_matrix_alloc(3, 3);
    gsl_matrix *V = gsl_matrix_alloc(3, 3);
    gsl_vector *S = gsl_vector_alloc(3);
    gsl_vector *work = gsl_vector_alloc(3);

    gsl_matrix_memcpy(U, H);
    gsl_linalg_SV_decomp(U, V, S, work);

    gsl_matrix_transpose(U);
    gsl_matrix *dM = gsl_matrix_alloc(3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, U, 0.0, dM);

    gsl_permutation *p = gsl_permutation_alloc(3);
    int signum;
    gsl_linalg_LU_decomp(dM, p, &signum);
    double d = sign(gsl_linalg_LU_det(dM, signum));

    gsl_matrix *Id = gsl_matrix_alloc(3, 3);
    gsl_matrix_set_identity(Id);
    gsl_matrix_set(Id, 2, 2, d);

    gsl_matrix *R = gsl_matrix_alloc(3, 3);
    gsl_matrix *R_temp = gsl_matrix_alloc(3, 3);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, Id, 0.0, R_temp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, R_temp, U, 0.0, R);

    return R;
}

double sign(double val)
{
    return (0.0 < val) - (val < 0.0);
}