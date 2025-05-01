//---------------------------------------------------------------------------
// discover.c
//
// Main SUBDUE discovery functions.
//
// SUBDUE 5
//---------------------------------------------------------------------------

#include "subdue.h"


//---------------------------------------------------------------------------
// NAME: DiscoverSubs
//
// INPUTS: (Parameters *parameters)
//
// RETURN: (SubList) - list of best discovered substructures
//
// PURPOSE: Discover the best substructures in the graphs according to
// the given parameters.  Note that we do not allow a single-vertex
// substructure of the form "SUB_#" on to the discovery list to avoid
// continually replacing "SUB_<n>" with "SUB_<n+1>".
//---------------------------------------------------------------------------

SubList *DiscoverSubs(Parameters *parameters)
{
   SubList *parentSubList;
   SubList *childSubList;
   SubList *extendedSubList;
   SubList *discoveredSubList;
   SubListNode *parentSubListNode;
   SubListNode *extendedSubListNode;
   Substructure *parentSub;
   Substructure *extendedSub;
   Substructure *recursiveSub = NULL;

   // parameters used
   ULONG limit          = parameters->limit;
   ULONG numBestSubs    = parameters->numBestSubs;
   ULONG beamWidth      = parameters->beamWidth;
   BOOLEAN valueBased   = parameters->valueBased;
   LabelList *labelList = parameters->labelList;
   BOOLEAN prune        = parameters->prune;
   ULONG maxVertices    = parameters->maxVertices;
   ULONG minVertices    = parameters->minVertices;
   ULONG outputLevel    = parameters->outputLevel;
   BOOLEAN recursion    = parameters->recursion;
   ULONG evalMethod     = parameters->evalMethod;

   // get initial one-vertex substructures
   // PrintGraph(parameters->posGraph, parameters->labelList);
   parentSubList = GetInitialSubs(parameters);

   discoveredSubList = AllocateSubList();
   while ((limit > 0) && (parentSubList->head != NULL)) 
   {
      parentSubListNode = parentSubList->head;
      childSubList = AllocateSubList();
      // extend each substructure in parent list
      while (parentSubListNode != NULL)
      {
         parentSub = parentSubListNode->sub;
         parentSubListNode->sub = NULL;
         if (outputLevel > 4) 
         {
            parameters->outputLevel = 1; // turn off instance printing
            printf("\nConsidering ");
            PrintSub(parentSub, parameters);
            printf("\n");
            parameters->outputLevel = outputLevel;
         }
         if ((((parentSub->numInstances > 1) && (evalMethod != EVAL_SETCOVER)) ||
              (parentSub->numNegInstances > 0)) &&
             (limit > 0))
         {
            limit--;
            if (outputLevel > 3)
               printf("%lu substructures left to be considered\n", limit);
            fflush(stdout);
            // InterconnectInstances(parentSub->instances, parameters->posGraph);
            // parentSub->definition = InstanceToGraph(parentSub->instances->head->instance, parameters->posGraph);
            // printf("Parent sub positions: ");
            // for(int vi=0; vi<parentSub->definition->numVertices; ++vi)
            // {
            //    printVector(parentSub->definition->vertices[vi].position,3);
            // }
            // printf("\n");
            extendedSubList = ExtendSub(parentSub, parameters);
            // printf("ds1\n");
            extendedSubListNode = extendedSubList->head;
            while (extendedSubListNode != NULL) 
            {
               // printf("Extended sub positions: ");
               // for(int vi=0; vi<extendedSubListNode->sub->definition->numVertices; ++vi)
               // {
               //    printVector(extendedSubListNode->sub->definition->vertices[vi].position,3);
               // }
               // printf("\n");
               extendedSub = extendedSubListNode->sub;
               // printf("ds1.1\n");
               extendedSubListNode->sub = NULL;
               // printf("Sub: num verts: %ld, max verts: %ld\n", extendedSub->definition->numVertices, maxVertices);
               if (extendedSub->definition->numVertices <= maxVertices) 
               {
                  // evaluate each extension and add to child list
                  // printf("ds1.2\n");
                  EvaluateSub(extendedSub, parameters);
                  // printf("ds1.2.0\n");
                  // PrintSub(extendedSub, parameters);
                  if (prune && (extendedSub->value < parentSub->value)) 
                  {
                     // printf("ds1.2.1\n");
                     FreeSub(extendedSub);
                  } 
                  else 
                  {
                     // printf("ds1.2.2\n");
                     SubListInsert(extendedSub, childSubList, beamWidth, 
                                   valueBased, labelList, parameters->use_pose_cost);
                  }
                  // printf("ds1.3\n");
               } 
               else 
               {
                  FreeSub(extendedSub);
               }
               extendedSubListNode = extendedSubListNode->next;
               // printf("ds1.4\n");
            }
            FreeSubList(extendedSubList);
         }
         // add parent substructure to final discovered list
         if (parentSub->definition->numVertices >= minVertices) 
         {
            if (! SinglePreviousSub(parentSub, parameters)) 
            {
               // consider recursive substructure, if requested
               if (recursion)
                  recursiveSub = RecursifySub(parentSub, parameters);
               if (outputLevel > 3)
                  PrintNewBestSub(parentSub, discoveredSubList, parameters);
               SubListInsert(parentSub, discoveredSubList, numBestSubs, FALSE,
                             labelList, parameters->use_pose_cost);
               if (recursion && (recursiveSub != NULL)) 
               {
                  if (outputLevel > 4) 
                  {
                     parameters->outputLevel = 1; // turn off instance printing
                     printf("\nConsidering Recursive ");
                     PrintSub(recursiveSub, parameters);
                     printf ("\n");
                     parameters->outputLevel = outputLevel;
                  }
                  if (outputLevel > 3)
                     PrintNewBestSub(recursiveSub, discoveredSubList, parameters);
                  SubListInsert(recursiveSub, discoveredSubList, numBestSubs,
                                FALSE, labelList, parameters->use_pose_cost);
               }
            }
         } 
         else 
         {
            FreeSub (parentSub);
         }
         parentSubListNode = parentSubListNode->next;
      }
      FreeSubList(parentSubList);
      parentSubList = childSubList;
      // printf("Best sub cost: %g\n", discoveredSubList->head->sub->instances->head->instance->minMatchCost);
   }

   if ((limit > 0) && (outputLevel > 2))
      printf ("\nSubstructure queue empty.\n");

   // try to insert any remaining subs in parent list on to discovered list
   parentSubListNode = parentSubList->head;
   while (parentSubListNode != NULL) 
   {
      parentSub = parentSubListNode->sub;
      parentSubListNode->sub = NULL;
      if (parentSub->definition->numVertices >= minVertices) 
      {
         if (! SinglePreviousSub(parentSub, parameters)) 
         {
            if (outputLevel > 3)
               PrintNewBestSub(parentSub, discoveredSubList, parameters);
            SubListInsert(parentSub, discoveredSubList, numBestSubs, FALSE,
                          labelList, parameters->use_pose_cost);
         }
      } 
      else 
      {
         FreeSub(parentSub);
      }
      parentSubListNode = parentSubListNode->next;
   }
   FreeSubList(parentSubList);
   // printf("Best sub cost: %g\n", discoveredSubList->head->sub->instances->head->instance->minMatchCost);
   return discoveredSubList;
}


//---------------------------------------------------------------------------
// NAME: GetInitialSubs
//
// INPUTS: (Parameters *parameters)
//
// RETURN: (SubList *)
//
// PURPOSE: Return a list of substructures, one for each unique vertex
// label in the positive graph that has at least two instances.
//---------------------------------------------------------------------------

SubList *GetInitialSubs(Parameters *parameters)
{
   SubList *initialSubs;
   ULONG i, j;
   ULONG vertexLabelIndex;
   ULONG numInitialSubs;
   Graph *g;
   Substructure *sub;
   Instance *instance;

   // parameters used
   Graph *posGraph      = parameters->posGraph;
   Graph *negGraph      = parameters->negGraph;
   LabelList *labelList = parameters->labelList;
   ULONG outputLevel    = parameters->outputLevel;
   ULONG currentIncrement = 0;
   ULONG startVertexIndex;

   if (parameters->incremental)
   {
      currentIncrement = GetCurrentIncrementNum(parameters);
      // Index for first vertex in increment
      // Begin with the index for the first vertex in this increment and
      // move up through all remaining vertices.  Relies on the fact that
      // each new increment is placed on the end of the vertex array and that
      // we are only interested in the current (last) increment
      startVertexIndex = GetStartVertexIndex(currentIncrement, parameters, POS);
      if (parameters->outputLevel > 2)
         printf("Start vertex index = %lu\n", startVertexIndex);
   }
   else 
      startVertexIndex = 0;

   // reset labels' used flag
   for (i = 0; i < labelList->numLabels; i++)
      labelList->labels[i].used = FALSE;
  
   numInitialSubs = 0;
   initialSubs = AllocateSubList();
   for (i = startVertexIndex; i < posGraph->numVertices; i++)
   {
      vertexLabelIndex = posGraph->vertices[i].label;
      if (labelList->labels[vertexLabelIndex].used == FALSE) 
      // if (labelList->labels[vertexLabelIndex].used == FALSE && labelList->labels[vertexLabelIndex].imp == TRUE) 
      {
         labelList->labels[vertexLabelIndex].used = TRUE;

         // create one-vertex substructure definition
         g = AllocateGraph(1, 0);
         g->vertices[0].label = vertexLabelIndex;
         g->vertices[0].numEdges = 0;
         g->vertices[0].edges = NULL;
         // g->vertices[0].position = gsl_vector_alloc(3);
         // g->vertices[0].orientation = gsl_vector_alloc(3);
         // gsl_vector_memcpy(g->vertices[0].position, posGraph->vertices[i].position);
         // gsl_vector_memcpy(g->vertices[0].orientation, posGraph->vertices[i].orientation);
         g->vertices[0].position = posGraph->vertices[i].position;
         g->vertices[0].orientation = posGraph->vertices[i].orientation;
         // allocate substructure
         sub = AllocateSub();
         sub->definition = g;
         sub->instances = AllocateInstanceList();
         // collect instances in positive graph
         j = posGraph->numVertices;
         do 
         {
            j--;
            if (posGraph->vertices[j].label == vertexLabelIndex) 
            {
               // ***** do inexact label matches here? (instance->minMatchCost
               // ***** too)
               instance = AllocateInstance(1, 0);
               instance->vertices[0] = j;
               instance->mapping[0].v1 = 0;
               instance->mapping[0].v2 = j;
               instance->minMatchCost = 0.0;
               InstanceListInsert(instance, sub->instances, FALSE);
               sub->numInstances++;
            }
         } while (j > i);

         // only keep substructure if more than one positive instance
         if (sub->numInstances > 1) 
         {
            if (negGraph != NULL) 
            {
               // collect instances in negative graph
               sub->negInstances = AllocateInstanceList();
               j = negGraph->numVertices;
               if (parameters->incremental)
                  startVertexIndex =
                     GetStartVertexIndex(currentIncrement, parameters, POS);
               else 
                  startVertexIndex = 0;
               do 
               {
                  j--;
                  if (negGraph->vertices[j].label == vertexLabelIndex) 
                  {
                     // ***** do inexact label matches here? 
                     // ***** (instance->minMatchCost too)
                     instance = AllocateInstance(1, 0);
                     instance->vertices[0] = j;
                     instance->mapping[0].v1 = 0;
                     instance->mapping[0].v2 = j;
                     instance->minMatchCost = 0.0;
                     InstanceListInsert(instance, sub->negInstances, FALSE);
                     sub->numNegInstances++;
                  }
                              // We need to try all negative graph labels
               } while (j > startVertexIndex);
            }
            EvaluateSub(sub, parameters);
            // add to initialSubs
            SubListInsert(sub, initialSubs, 0, FALSE, labelList, parameters->use_pose_cost);
            numInitialSubs++;
         } 
         else 
         { // prune single-instance substructure
            FreeSub(sub);
         }
      }
   }
   if (outputLevel > 1)
   {
      printf("%lu initial substructures\n", numInitialSubs);
      // PrintSubList(initialSubs, parameters);
   }


   return initialSubs;
}


//---------------------------------------------------------------------------
// NAME: SinglePreviousSub
//
// INPUTS: (Substructure *sub) - substructure to check
//         (Parameters *parameters)
//
// RETURN: (BOOLEAN)
//
// PURPOSE: Returns TRUE if the given substructure is a single-vertex
// substructure and the vertex refers to a previously-discovered
// substructure, i.e., the vertex label is of the form "SUB_#".  This
// is used to prevent repeatedly compressing the graph by replacing a
// "SUB_<n>" vertex by a "SUB_<n+1>" vertex.
//---------------------------------------------------------------------------

BOOLEAN SinglePreviousSub(Substructure *sub, Parameters *parameters)
{
   BOOLEAN match;
   // parameters used
   LabelList *labelList = parameters->labelList;

   match = FALSE;
   if ((sub->definition->numVertices == 1) &&  // single-vertex sub?
       (SubLabelNumber (sub->definition->vertices[0].label, labelList) > 0))
       // valid substructure label
      match = TRUE;

   return match;
}
