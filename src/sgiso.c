//---------------------------------------------------------------------------
// sgiso.c
//
// These functions implement a subgraph isomorphism algorithm embodied
// in the FindInstances function.  This is primarily used to find
// predefined substructures.
//
// Subdue 5
//---------------------------------------------------------------------------

#include "subdue.h"


//---------------------------------------------------------------------------
// NAME: FindInstances
//
// INPUTS: (Graph *g1) - graph to search for
//         (Graph *g2) - graph to search in
//         (Parameters *parameters)
//
// RETURN: (InstanceList *) - list of instances of g1 in g2, may be empty
//
// PURPOSE: Searches for subgraphs of g2 that match g1 and returns the
// list of such subgraphs as instances in g2.  Returns empty list if
// no matches exist.  This procedure mimics the DiscoverSubs loop by
// repeatedly expanding instances of subgraphs of g1 in g2 until
// matches are found.  The procedure is optimized toward g1 being a
// small graph and g2 being a large graph.
//
// Note: This procedure is equivalent to the NP-Hard subgraph
// isomorphism problem, and therefore can be quite slow for some
// graphs.
//---------------------------------------------------------------------------

InstanceList *FindInstances(Graph *g1, Graph *g2, Parameters *parameters)
{
   ULONG v;
   ULONG v1;
   ULONG e1;
   Vertex *vertex1;
   Edge *edge1;
   InstanceList *instanceList;
   BOOLEAN *reached;
   BOOLEAN noMatches;
   BOOLEAN found;

   reached = (BOOLEAN *) malloc(sizeof(BOOLEAN) * g1->numVertices);
   if (reached == NULL)
      OutOfMemoryError("FindInstances:reached");
   for (v1 = 0; v1 < g1->numVertices; v1++)
      reached[v1] = FALSE;
   v1 = 0; // first vertex in g1
   reached[v1] = TRUE;
   vertex1 = & g1->vertices[v1];
   instanceList = FindSingleVertexInstances(g2, vertex1, parameters);
   noMatches = FALSE;
   if (instanceList->head == NULL) // no matches to vertex1 found in g2
      noMatches = TRUE;
   while ((vertex1 != NULL) && (! noMatches)) 
   {
      vertex1->used = TRUE;
      // extend by each unmarked edge involving vertex v1
      for (e1 = 0; ((e1 < vertex1->numEdges) && (! noMatches)); e1++) 
      {
         edge1 = & g1->edges[g1->vertices[v1].edges[e1]];
         if (! edge1->used) 
         {
            reached[edge1->vertex1] = TRUE;
            reached[edge1->vertex2] = TRUE;
            edge1->used = TRUE;
            instanceList =
               ExtendInstancesByEdge(instanceList, g1, edge1, g2, parameters);
            if (instanceList->head == NULL)
               noMatches = TRUE;
         }
      }
      // find next un-used, reached vertex
      vertex1 = NULL;
      found = FALSE;
      for (v = 0; ((v < g1->numVertices) && (! found)); v++) 
      {
         if ((! g1->vertices[v].used) && (reached[v])) 
         {
            v1 = v;
            vertex1 = & g1->vertices[v1];
            found = TRUE;
         }
      }
   }
   free(reached);

   // set used flags to FALSE
   for (v1 = 0; v1 < g1->numVertices; v1++)
      g1->vertices[v1].used = FALSE;
   for (e1 = 0; e1 < g1->numEdges; e1++)
      g1->edges[e1].used = FALSE;

   // filter instances not matching g1
   // filter overlapping instances if appropriate
   instanceList = FilterInstances(g1, instanceList, g2, parameters);
 
   return instanceList;
}


Instance *FindFittingInstances(Graph *g1, Graph *g2, Parameters *parameters)
{
   ULONG v;
   ULONG v1;
   ULONG e1;
   Vertex *vertex1;
   Edge *edge1;
   InstanceList *instanceList;
   InstanceListNode * instance_list_node;
   Instance * best_matching_instance;
   double best_fitness_score = 0.0;
   double fitness_thr = 0.7;
   BOOLEAN *reached;
   BOOLEAN noMatches;
   BOOLEAN found;
   ULONG pattern_size;
   ULONG subgraph_size;
   pattern_size = g1->numVertices + g1->numEdges;

   reached = (BOOLEAN *) malloc(sizeof(BOOLEAN) * g1->numVertices);
   if (reached == NULL)
      OutOfMemoryError("FindInstances:reached");
   for (v1 = 0; v1 < g1->numVertices; v1++)
      reached[v1] = FALSE;
   v1 = 0; // !!!! Change this to use a vertex that has the highest probability of being in the pattern
   reached[v1] = TRUE;
   vertex1 = & g1->vertices[v1];
   instanceList = FindSingleVertexInstances(g2, vertex1, parameters);
   noMatches = FALSE;
   if (instanceList->head == NULL) // no matches to vertex1 found in g2 --> In this case try other vertex of same class
      noMatches = TRUE;
   while ((vertex1 != NULL) && (! noMatches)) 
   {
      vertex1->used = TRUE;
      // extend by each unmarked edge involving vertex v1
      for (e1 = 0; ((e1 < vertex1->numEdges) && (! noMatches)); e1++) 
      {
         edge1 = & g1->edges[g1->vertices[v1].edges[e1]];
         // printf("------\n");
         // printf("Trying for edge %s - %s | %s\n", 
         //    parameters->labelList->labels[g1->vertices[edge1->vertex1].label].labelValue.stringLabel, 
         //    parameters->labelList->labels[g1->vertices[edge1->vertex2].label].labelValue.stringLabel,
         //    parameters->labelList->labels[edge1->label].labelValue.stringLabel);
         if (! edge1->used) 
         {
            reached[edge1->vertex1] = TRUE;
            reached[edge1->vertex2] = TRUE;
            edge1->used = TRUE;
            instanceList =
               ExtendInstancesByEdgeSoft(instanceList, g1, edge1, g2, parameters, reached);
            if (instanceList->head == NULL)
               noMatches = TRUE;

            if(instanceList->head != NULL)
            {
               // printf("Instance scores: ");
               instance_list_node = instanceList->head;
               while(instance_list_node != NULL)
               {
                  double fitness_score = (double)(instance_list_node->instance->numVertices + instance_list_node->instance->numEdges) / pattern_size;
                  // printf("%g, ", fitness_score);
                  if(fitness_score > best_fitness_score)
                  {
                     best_fitness_score = fitness_score;
                     best_matching_instance = CopyInstance(instance_list_node->instance);
                  }
                  instance_list_node = instance_list_node->next;
               }
               // printf("\n");
            }
         }
      }
      // printf("Best fitting score so far: %g\n", best_fitness_score);
      
      if(best_fitness_score >= fitness_thr)
      {
         break;
      }
            
      // find next un-used, reached vertex
      vertex1 = NULL;
      found = FALSE;
      for (v = 0; ((v < g1->numVertices) && (! found)); v++) 
      {
         if ((! g1->vertices[v].used) && (reached[v])) 
         {
            v1 = v;
            vertex1 = & g1->vertices[v1];
            found = TRUE;
         }
      }
   }
   free(reached);

   // set used flags to FALSE
   for (v1 = 0; v1 < g1->numVertices; v1++)
      g1->vertices[v1].used = FALSE;
   for (e1 = 0; e1 < g1->numEdges; e1++)
      g1->edges[e1].used = FALSE;

   // filter instances not matching g1
   // filter overlapping instances if appropriate
   instanceList = FilterInstances(g1, instanceList, g2, parameters);
 
   return best_matching_instance;
}



//---------------------------------------------------------------------------
// NAME: FindSingleVertexInstances
//
// INPUTS: (Graph *graph) - graph to look in for matching vertices
//         (Vertex *vertex) - vertex to look for
//         (Parameters *parameters)
//
// RETURN: (InstanceList *) - list of single-vertex instances matching
//                            vertex
//
// PURPOSE: Return a (possibly-empty) list of single-vertex instances, one
// for each vertex in the given graph that matches the given vertex.
//---------------------------------------------------------------------------

InstanceList *FindSingleVertexInstances(Graph *graph, Vertex *vertex,
                                        Parameters *parameters)
{
   ULONG v;
   InstanceList *instanceList;
   Instance *instance;

   instanceList = AllocateInstanceList();
   for (v = 0; v < graph->numVertices; v++) 
   {
      if (graph->vertices[v].label == vertex->label && !graph->vertices[v].in_pattern) 
      {
         // ***** do inexact label matches here? (instance->minMatchCost too)
         instance = AllocateInstance(1, 0);
         instance->vertices[0] = v;
         instance->minMatchCost = 0.0;
         InstanceListInsert(instance, instanceList, FALSE);
      }
   }
   return instanceList;
}


//---------------------------------------------------------------------------
// NAME: ExtendInstancesByEdge
//
// INPUTS: (InstanceList *instanceList) - instances to extend by one edge
//         (Graph *g1) - graph whose instances we are looking for
//         (Edge *edge1) - edge in g1 by which to extend each instance
//         (Graph *g2) - graph containing instances
//         (Parameters parameters)
//
// RETURN: (InstanceList *) - new instance list with extended instances
//
// PURPOSE: Attempts to extend each instance in instanceList by an
// edge from graph g2 that matches the attributes of the given edge in
// graph g1.  Returns a new (possibly empty) instance list containing
// the extended instances.  The given instance list is de-allocated.
//---------------------------------------------------------------------------

InstanceList *ExtendInstancesByEdge(InstanceList *instanceList,
                                    Graph *g1, Edge *edge1, Graph *g2,
                                    Parameters *parameters)
{
   InstanceList *newInstanceList;
   InstanceListNode *instanceListNode;
   Instance *instance;
   Instance *newInstance;
   ULONG v2;
   ULONG e2;
   Edge *edge2;
   Vertex *vertex2;

   newInstanceList = AllocateInstanceList();
   // extend each instance
   instanceListNode = instanceList->head;
   while (instanceListNode != NULL) 
   {
      instance = instanceListNode->instance;
      MarkInstanceEdges(instance, g2, TRUE);
      // consider extending from each vertex in instance
      for (v2 = 0; v2 < instance->numVertices; v2++) 
      {
         vertex2 = & g2->vertices[instance->vertices[v2]];
         for (e2 = 0; e2 < vertex2->numEdges; e2++) 
         {
            edge2 = & g2->edges[vertex2->edges[e2]];
            if ((! edge2->used) &&
                (EdgesMatch(g1, edge1, g2, edge2, parameters))) 
            {
               // add new instance to list
               newInstance =
                  CreateExtendedInstance(instance, instance->vertices[v2],
                                         vertex2->edges[e2], g2);
               InstanceListInsert(newInstance, newInstanceList, TRUE);
            }
         }
      }
      MarkInstanceEdges(instance, g2, FALSE);
      instanceListNode = instanceListNode->next;
   }
   FreeInstanceList(instanceList);
   return newInstanceList;
}

//---------------------------------------------------------------------------
// NAME: ExtendInstancesByEdgeSoft
//
// INPUTS: (InstanceList *instanceList) - instances to extend by one edge
//         (Graph *g1) - graph whose instances we are looking for
//         (Edge *edge1) - edge in g1 by which to extend each instance
//         (Graph *g2) - graph containing instances
//         (Parameters parameters)
//
// RETURN: (InstanceList *) - new instance list with extended instances
//
// PURPOSE: Does the same thing as ExtendInstancesByEdge but will append the 
// instances that were not expanded by this edge to the output list.
//---------------------------------------------------------------------------
InstanceList *ExtendInstancesByEdgeSoft(InstanceList *instanceList,
                                    Graph *g1, Edge *edge1, Graph *g2,
                                    Parameters *parameters, BOOLEAN *reached)
{
   InstanceList *newInstanceList;
   InstanceListNode *instanceListNode;
   Instance *instance;
   Instance *newInstance;
   ULONG v2;
   ULONG e2;
   Edge *edge2;
   Vertex *vertex2;

   newInstanceList = AllocateInstanceList();
   // extend each instance
   instanceListNode = instanceList->head;
   int num_instances = 0;
   while (instanceListNode != NULL) 
   {
      ++num_instances;
      // printf("Instance %d\n", num_instances);
      instance = instanceListNode->instance;
      MarkInstanceEdges(instance, g2, TRUE);
      // consider extending from each vertex in instance
      BOOLEAN found_one_extension = FALSE;
      BOOLEAN reached_max_expansion = TRUE;
      // printf("Num vertices: %lu \n", instance->numVertices);
      // PrintInstance(instance, 0, g2, parameters->labelList);
      for (v2 = 0; v2 < instance->numVertices; v2++) 
      {
         vertex2 = & g2->vertices[instance->vertices[v2]];
         // printf("Vertex %lu: Num edges: %lu \n", v2, vertex2->numEdges);
         for (e2 = 0; e2 < vertex2->numEdges; e2++) 
         {
            edge2 = & g2->edges[vertex2->edges[e2]];
            if ((! edge2->used))
            {
               reached_max_expansion = FALSE;
               if((EdgesMatch(g1, edge1, g2, edge2, parameters))) 
               {
                  ULONG other_vertex_ind;
                  if(edge2->vertex1 == instance->vertices[v2])
                     other_vertex_ind = edge2->vertex2;
                  else
                     other_vertex_ind = edge2->vertex1;
                  
                  if(!g2->vertices[other_vertex_ind].in_pattern)
                  {
                     // add new instance to list
                     // printf("Instance before extending:\n");
                     // PrintInstance(instance, 0, g2, parameters->labelList);
                     newInstance =
                        CreateExtendedInstance(instance, instance->vertices[v2],
                                             vertex2->edges[e2], g2);
                     if(InstanceSubgraphMatch(newInstance, g1, g2, parameters, reached))
                     {
                        InstanceListInsert(newInstance, newInstanceList, TRUE);
                        found_one_extension = TRUE;
                     }   
                  }
               }
            }
         }
      }
      if(!found_one_extension && !reached_max_expansion)
      {
         // printf("Could not extend instance\n");
         newInstance = CopyInstance(instance);
         InstanceListInsert(newInstance, newInstanceList, TRUE);
      }
      MarkInstanceEdges(instance, g2, FALSE);
      instanceListNode = instanceListNode->next;
   }
   // printf("Num instances expanded: %d\n", num_instances);
   FreeInstanceList(instanceList);
   return newInstanceList;
}


//---------------------------------------------------------------------------
// NAME: EdgesMatch
//
// INPUTS: (Graph *g1) - graph containing edge e1
//         (Edge *edge1) - edge from g1 to be matched
//         (Graph *g2) - graph containing edge e2
//         (Edge *edge2) - edge from g2 being matched
//         (Parameters *paramters)
//
// RETURN: (BOOLEAN) - TRUE if edges and their vertices match; else FALSE
//
// PURPOSE: Check that the given edges match in terms of their labels,
// directedness, and their starting and ending vertex labels.
//
// ***** do inexact label matches?
//---------------------------------------------------------------------------

BOOLEAN EdgesMatch(Graph *g1, Edge *edge1, Graph *g2, Edge *edge2,
                   Parameters *parameters)
{
   ULONG vertex11Label;
   ULONG vertex12Label;
   ULONG vertex21Label;
   ULONG vertex22Label;
   BOOLEAN match;

   match = FALSE;
   vertex11Label = g1->vertices[edge1->vertex1].label;
   vertex12Label = g1->vertices[edge1->vertex2].label;
   vertex21Label = g2->vertices[edge2->vertex1].label;
   vertex22Label = g2->vertices[edge2->vertex2].label;
   if ((edge1->label == edge2->label) &&
       (edge1->directed == edge2->directed)) 
   {
      if ((vertex11Label == vertex21Label) &&
          (vertex12Label == vertex22Label))
         match = TRUE;
      else if ((! edge1->directed) &&
               (vertex11Label == vertex22Label) &&
               (vertex12Label == vertex21Label))
         match = TRUE;
   }
   return match;
}


//---------------------------------------------------------------------------
// NAME: FilterInstances
//
// INPUTS: (Graph *subGraph) - graph that instances must match
//         (InstanceList *instanceList) - list of instances
//         (Graph *graph) - graph containing instances
//         (Parameters *parameters)
//
// RETURN: (InstanceList *) - filtered instance list
//
// PURPOSE: Creates and returns a new instance list containing only
// those instances matching subGraph.  If
// parameters->allowInstanceOverlap=FALSE, then remaining instances
// will not overlap.  The given instance list is de-allocated.
//---------------------------------------------------------------------------

InstanceList *FilterInstances(Graph *subGraph, InstanceList *instanceList,
                              Graph *graph, Parameters *parameters)
{
   InstanceListNode *instanceListNode;
   InstanceListNode *overlapping_instances;
   Instance *instance;
   InstanceList *newInstanceList;
   Graph *instanceGraph;
   double thresholdLimit;
   double matchCost;

   newInstanceList = AllocateInstanceList();
   if (instanceList != NULL) 
   {
      instanceListNode = instanceList->head;
      while (instanceListNode != NULL) 
      {
         if (instanceListNode->instance != NULL) 
         {
            instance = instanceListNode->instance;
            int overlap_index;
            if (parameters->allowInstanceOverlap ||
                (! InstanceListOverlap(instance, newInstanceList, &overlapping_instances, &overlap_index))) 
            {
               thresholdLimit = parameters->threshold *
                                (instance->numVertices + instance->numEdges);
               instanceGraph = InstanceToGraph(instance, graph);
               if (GraphMatch(subGraph, instanceGraph, parameters->labelList,
                              thresholdLimit, & matchCost, NULL, parameters->use_pose_cost)) 
               {
                  if (matchCost < instance->minMatchCost)
                     instance->minMatchCost = matchCost;
                  InstanceListInsert(instance, newInstanceList, FALSE);
               }
               FreeGraph(instanceGraph);
            }
         }
         instanceListNode = instanceListNode->next;
      }
   }
   FreeInstanceList(instanceList);
   FreeInstanceListNode(overlapping_instances);
   return newInstanceList;
}


BOOLEAN InstanceSubgraphMatch(Instance *instance, Graph *subgraph, Graph *g2, Parameters *parameters, BOOLEAN *reached)
{
   double thresholdLimit = parameters->threshold * (instance->numVertices + instance->numEdges);
   Graph *instance_graph;
   Graph *subset_graph;
   double match_cost;
   BOOLEAN match = FALSE;
   Edge *edge;
   Vertex *vertex;
   ULONG *vmap;
   int map_size = 0;

   vmap = (ULONG *)malloc(sizeof(ULONG) * subgraph->numVertices);

   subset_graph = AllocateGraph(0,0);
   for(ULONG v=0; v<subgraph->numVertices; ++v)
   {
      vertex = & subgraph->vertices[v];
      if(reached[v])
      {
         AddVertex(subset_graph, vertex->label);
         subset_graph->vertices[subset_graph->numVertices-1].position = vertex->position;
         vmap[v] = subset_graph->numVertices - 1;
         ++map_size;
      }
   }
   for(ULONG e=0; e<subgraph->numEdges; ++e)
   {
      edge = & subgraph->edges[e];
      // if(reached[edge->vertex1] && reached[edge->vertex2])
      if(edge->used)
      {
         AddEdge(subset_graph, vmap[edge->vertex1], vmap[edge->vertex2], edge->directed, edge->label, FALSE);
         edge->used = TRUE;
      }
   }

   // printf("vmap: ");
   // for(int i=0; i<subgraph->numVertices; ++i)
   // {
   //    if(reached[i])
   //    {
   //       printf("%lu : %lu | ", i, vmap[i]);
   //    }
   // }
   // printf("\n");

   // printf("Subgraph:\n");
   // PrintGraph(subgraph, parameters->labelList);
   // printf("Subset graph:\n");
   // PrintGraph(subset_graph, parameters->labelList);

   instance_graph = InstanceToGraph(instance, g2);
   // printf("instance graph:\n");
   // PrintGraph(instance_graph, parameters->labelList);
   if (GraphMatch(subset_graph, instance_graph, parameters->labelList,
                  thresholdLimit, & match_cost, NULL, parameters->use_pose_cost)) 
   {
      match = TRUE;
   }
   // printf("Match: %d\n", match);
   
   FreeGraph(instance_graph);

   return match;
}
