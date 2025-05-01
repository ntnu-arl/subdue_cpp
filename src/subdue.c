#include "subdue.h"

#include "time.h"
#include <sys/times.h>
#include <unistd.h>

Parameters *initializeNewParameters()
{
   Parameters *params = (Parameters *)malloc(sizeof(Parameters));
   if (params == NULL)
      OutOfMemoryError("parameters");
   return params;
}

Parameters *processParameters(Parameters *parameters)
{
   parameters->log2Factorial = (double *)malloc(2 * sizeof(double));
   if (parameters->log2Factorial == NULL)
      OutOfMemoryError("GetParameters:parameters->log2Factorial");
   parameters->log2FactorialSize = 2;
   parameters->log2Factorial[0] = 0; // lg(0!)
   parameters->log2Factorial[1] = 0; // lg(1!)

   // read graphs from input file
   parameters->labelList = AllocateLabelList();
   parameters->posGraph = NULL;
   parameters->negGraph = NULL;
   parameters->numPosEgs = 0;
   parameters->numNegEgs = 0;
   parameters->posEgsVertexIndices = NULL;
   parameters->negEgsVertexIndices = NULL;

   if (parameters->incremental)
   {
      if (parameters->predefinedSubs)
      {
         fprintf(stderr, "Cannot process predefined examples incrementally");
         exit(1);
      }

      if (parameters->evalMethod == EVAL_MDL)
      {
         fprintf(stderr, "Incremental SUBDUE does not support EVAL_MDL, ");
         fprintf(stderr, "switching to EVAL_SIZE\n");
         parameters->evalMethod = EVAL_SIZE;
      }

      if ((parameters->evalMethod == EVAL_SIZE) && (parameters->compress))
      {
         fprintf(stderr, "Incremental SUBDUE does not support compression, ");
         fprintf(stderr, "with EVAL_SIZE, turning compression off\n");
         parameters->compress = FALSE;
      }

      if (parameters->iterations > 1)
      {
         fprintf(stderr,
                 "Incremental SUBDUE only one iteration, setting to 1\n");
         parameters->iterations = 1;
      }
   }

   // read predefined substructures
   parameters->numPreSubs = 0;
   if (parameters->predefinedSubs)
      ReadPredefinedSubsFile(parameters);

   parameters->incrementList = (IncrementList *)malloc(sizeof(IncrementList));
   parameters->incrementList->head = NULL;

   if (parameters->incremental)
   {
      parameters->vertexList = (InstanceVertexList *)malloc(sizeof(InstanceVertexList));
      parameters->vertexList->avlTreeList = (AvlTreeList *)malloc(sizeof(AvlTreeList));
      parameters->vertexList->avlTreeList->head = NULL;
   }

   // create output file, if given
   FILE *outputFile;
   if (parameters->outputToFile)
   {
      outputFile = fopen(parameters->outFileName, "w");
      if (outputFile == NULL)
      {
         printf("ERROR: unable to write to output file %s\n",
                parameters->outFileName);
         exit(1);
      }
      fclose(outputFile);
   }

   /// PostProcessParameters
   Increment *increment = NULL;

   if (parameters->incremental)
      increment = GetCurrentIncrement(parameters);

   // Code from this point until end of function was moved from GetParameters
   // if (parameters->numPosEgs == 0)
   // {
   //    fprintf(stderr, "ERROR: no positive graphs defined\n");
   //    exit(1);
   // }

   // Check bounds on discovered substructures' number of vertices
   // if (parameters->maxVertices == 0)
   //    parameters->maxVertices = parameters->posGraph->numVertices;
   // if (parameters->maxVertices < parameters->minVertices)
   // {
   //    fprintf(stderr, "ERROR: minsize exceeds maxsize\n");
   //    exit(1);
   // }

   // Set limit accordingly
   if (parameters->limit == 0)
   {
      if (parameters->incremental)
         parameters->limit = increment->numPosEdges / 2;
      else
         parameters->limit = parameters->posGraph->numEdges / 2;
   }

   return parameters;
}

//---------------------------------------------------------------------------
// NAME: GetParameters
//
// INPUTS: (int argc) - number of command-line arguments
//         (char *argv[]) - array of command-line argument strings
//
// RETURN: (Parameters *)
//
// PURPOSE: Initialize parameters structure and process command-line
//          options.
//---------------------------------------------------------------------------

Parameters *GetParameters(int argc, char **argv)
{
   Parameters *parameters;
   int i;
   double doubleArg;
   ULONG ulongArg;
   BOOLEAN limitSet = FALSE;
   FILE *outputFile;

   parameters = (Parameters *)malloc(sizeof(Parameters));
   if (parameters == NULL)
      OutOfMemoryError("parameters");

   // initialize default parameter settings
   parameters->directed = TRUE;
   parameters->limit = 0;
   parameters->numBestSubs = 3;
   parameters->beamWidth = 4;
   parameters->valueBased = FALSE;
   parameters->prune = FALSE;
   strcpy(parameters->outFileName, "none");
   parameters->outputToFile = FALSE;
   parameters->outputLevel = 2;
   parameters->allowInstanceOverlap = FALSE;
   parameters->threshold = 0.0;
   parameters->evalMethod = EVAL_MDL;
   parameters->iterations = 1;
   strcpy(parameters->psInputFileName, "none");
   parameters->predefinedSubs = FALSE;
   parameters->minVertices = 1;
   parameters->maxVertices = 0; // i.e., infinity
   parameters->recursion = FALSE;
   parameters->variables = FALSE;
   parameters->relations = FALSE;
   parameters->incremental = FALSE;
   parameters->compress = FALSE;

   if (argc < 2)
   {
      fprintf(stderr, "input graph file name must be supplied\n");
      exit(1);
   }

   // process command-line options
   i = 1;
   while (i < (argc - 1))
   {
      if (strcmp(argv[i], "-beam") == 0)
      {
         i++;
         sscanf(argv[i], "%lu", &ulongArg);
         if (ulongArg == 0)
         {
            fprintf(stderr, "%s: beam must be greater than zero\n", argv[0]);
            exit(1);
         }
         parameters->beamWidth = ulongArg;
      }
      else if (strcmp(argv[i], "-compress") == 0)
      {
         parameters->compress = TRUE;
      }
      else if (strcmp(argv[i], "-eval") == 0)
      {
         i++;
         sscanf(argv[i], "%lu", &ulongArg);
         if ((ulongArg < 1) || (ulongArg > 3))
         {
            fprintf(stderr, "%s: eval must be 1-3\n", argv[0]);
            exit(1);
         }
         parameters->evalMethod = ulongArg;
      }
      else if (strcmp(argv[i], "-inc") == 0)
      {
         parameters->incremental = TRUE;
      }
      else if (strcmp(argv[i], "-iterations") == 0)
      {
         i++;
         sscanf(argv[i], "%lu", &ulongArg);
         parameters->iterations = ulongArg;
      }
      else if (strcmp(argv[i], "-limit") == 0)
      {
         i++;
         sscanf(argv[i], "%lu", &ulongArg);
         if (ulongArg == 0)
         {
            fprintf(stderr, "%s: limit must be greater than zero\n", argv[0]);
            exit(1);
         }
         parameters->limit = ulongArg;
         limitSet = TRUE;
      }
      else if (strcmp(argv[i], "-maxsize") == 0)
      {
         i++;
         sscanf(argv[i], "%lu", &ulongArg);
         if (ulongArg == 0)
         {
            fprintf(stderr, "%s: maxsize must be greater than zero\n", argv[0]);
            exit(1);
         }
         parameters->maxVertices = ulongArg;
      }
      else if (strcmp(argv[i], "-minsize") == 0)
      {
         i++;
         sscanf(argv[i], "%lu", &ulongArg);
         if (ulongArg == 0)
         {
            fprintf(stderr, "%s: minsize must be greater than zero\n", argv[0]);
            exit(1);
         }
         parameters->minVertices = ulongArg;
      }
      else if (strcmp(argv[i], "-nsubs") == 0)
      {
         i++;
         sscanf(argv[i], "%lu", &ulongArg);
         if (ulongArg == 0)
         {
            fprintf(stderr, "%s: nsubs must be greater than zero\n", argv[0]);
            exit(1);
         }
         parameters->numBestSubs = ulongArg;
      }
      else if (strcmp(argv[i], "-out") == 0)
      {
         i++;
         strcpy(parameters->outFileName, argv[i]);
         parameters->outputToFile = TRUE;
      }
      else if (strcmp(argv[i], "-output") == 0)
      {
         i++;
         sscanf(argv[i], "%lu", &ulongArg);
         if ((ulongArg < 1) || (ulongArg > 5))
         {
            fprintf(stderr, "%s: output must be 1-5\n", argv[0]);
            exit(1);
         }
         parameters->outputLevel = ulongArg;
      }
      else if (strcmp(argv[i], "-overlap") == 0)
      {
         parameters->allowInstanceOverlap = TRUE;
      }
      else if (strcmp(argv[i], "-prune") == 0)
      {
         parameters->prune = TRUE;
      }
      else if (strcmp(argv[i], "-ps") == 0)
      {
         i++;
         strcpy(parameters->psInputFileName, argv[i]);
         parameters->predefinedSubs = TRUE;
      }
      else if (strcmp(argv[i], "-recursion") == 0)
      {
         parameters->recursion = TRUE;
      }
      else if (strcmp(argv[i], "-relations") == 0)
      {
         parameters->relations = TRUE;
         parameters->variables = TRUE; // relations must involve variables
      }
      else if (strcmp(argv[i], "-threshold") == 0)
      {
         i++;
         sscanf(argv[i], "%lf", &doubleArg);
         if ((doubleArg < (double)0.0) || (doubleArg > (double)1.0))
         {
            fprintf(stderr, "%s: threshold must be 0.0-1.0\n", argv[0]);
            exit(1);
         }
         parameters->threshold = doubleArg;
      }
      else if (strcmp(argv[i], "-undirected") == 0)
      {
         parameters->directed = FALSE;
      }
      else if (strcmp(argv[i], "-valuebased") == 0)
      {
         parameters->valueBased = TRUE;
      }
      else if (strcmp(argv[i], "-variables") == 0)
      {
         parameters->variables = TRUE;
      }
      else
      {
         fprintf(stderr, "%s: unknown option %s\n", argv[0], argv[i]);
         exit(1);
      }
      i++;
   }

   if (parameters->iterations == 0)
      parameters->iterations = MAX_UNSIGNED_LONG; // infinity

   // initialize log2Factorial[0..1]
   parameters->log2Factorial = (double *)malloc(2 * sizeof(double));
   if (parameters->log2Factorial == NULL)
      OutOfMemoryError("GetParameters:parameters->log2Factorial");
   parameters->log2FactorialSize = 2;
   parameters->log2Factorial[0] = 0; // lg(0!)
   parameters->log2Factorial[1] = 0; // lg(1!)

   // read graphs from input file
   strcpy(parameters->inputFileName, argv[argc - 1]);
   parameters->labelList = AllocateLabelList();
   parameters->posGraph = NULL;
   parameters->negGraph = NULL;
   parameters->numPosEgs = 0;
   parameters->numNegEgs = 0;
   parameters->posEgsVertexIndices = NULL;
   parameters->negEgsVertexIndices = NULL;

   if (parameters->incremental)
   {
      if (parameters->predefinedSubs)
      {
         fprintf(stderr, "Cannot process predefined examples incrementally");
         exit(1);
      }

      if (parameters->evalMethod == EVAL_MDL)
      {
         fprintf(stderr, "Incremental SUBDUE does not support EVAL_MDL, ");
         fprintf(stderr, "switching to EVAL_SIZE\n");
         parameters->evalMethod = EVAL_SIZE;
      }

      if ((parameters->evalMethod == EVAL_SIZE) && (parameters->compress))
      {
         fprintf(stderr, "Incremental SUBDUE does not support compression, ");
         fprintf(stderr, "with EVAL_SIZE, turning compression off\n");
         parameters->compress = FALSE;
      }

      if (parameters->iterations > 1)
      {
         fprintf(stderr,
                 "Incremental SUBDUE only one iteration, setting to 1\n");
         parameters->iterations = 1;
      }
   }
   // else
   // {
   // ReadInputFile(parameters);
   // if (parameters->evalMethod == EVAL_MDL)
   // {
   //    parameters->posGraphDL = MDL(parameters->posGraph,
   //                                 parameters->labelList->numLabels, parameters);
   //    if (parameters->negGraph != NULL)
   //    {
   //       parameters->negGraphDL =
   //           MDL(parameters->negGraph, parameters->labelList->numLabels,
   //               parameters);
   //    }
   // }
   // }

   // read predefined substructures
   parameters->numPreSubs = 0;
   if (parameters->predefinedSubs)
      ReadPredefinedSubsFile(parameters);

   parameters->incrementList = malloc(sizeof(IncrementList));
   parameters->incrementList->head = NULL;

   if (parameters->incremental)
   {
      parameters->vertexList = malloc(sizeof(InstanceVertexList));
      parameters->vertexList->avlTreeList = malloc(sizeof(AvlTreeList));
      parameters->vertexList->avlTreeList->head = NULL;
   }

   // create output file, if given
   if (parameters->outputToFile)
   {
      outputFile = fopen(parameters->outFileName, "w");
      if (outputFile == NULL)
      {
         printf("ERROR: unable to write to output file %s\n",
                parameters->outFileName);
         exit(1);
      }
      fclose(outputFile);
   }

   return parameters;
}

void FreeParameters(Parameters *parameters)
{
   FreeGraph(parameters->posGraph);
   FreeGraph(parameters->negGraph);
   FreeLabelList(parameters->labelList);
   free(parameters->posEgsVertexIndices);
   free(parameters->negEgsVertexIndices);
   free(parameters->log2Factorial);
   free(parameters);
}

void runSubdue(Parameters *parameters, Graph *compressed_graphs[], LabelList *all_label_lists[], Substructure *discovered_subs[])
{
   struct tms tmsstart, tmsend;
   clock_t startTime, endTime;
   static long clktck = 0;
   time_t iterationStartTime;
   time_t iterationEndTime;
   SubList *subList;
   FILE *outputFile;
   ULONG iteration;
   BOOLEAN done;

   char sub_names[parameters->iterations][TOKEN_LEN];
   // Graph **discovered_subs;
   // LabelList **all_label_lists;

   // discovered_subs = (Graph **)realloc(discovered_subs, (sizeof(Graph *) * parameters->iterations));
   // all_label_lists = (LabelList **)realloc(all_label_lists, (sizeof(LabelList *) * parameters->iterations));
   // compressed_graphs = (Graph **)realloc(compressed_graphs, (sizeof(Graph *) * parameters->iterations));

   clktck = sysconf(_SC_CLK_TCK);
   startTime = times(&tmsstart);

   // compress pos and neg graphs with predefined subs, if given
   if (parameters->numPreSubs > 0)
      CompressWithPredefinedSubs(parameters);

   // PostProcessParameters(parameters);
   // PrintParameters(parameters);

   if (parameters->iterations > 1)
      printf("----- Iteration 1 -----\n\n");

   
   iteration = 1;
   done = FALSE;
   while ((iteration <= parameters->iterations) && (!done))
   {
      iterationStartTime = time(NULL);
      if (iteration > 1)
         printf("----- Iteration %lu -----\n\n", iteration);

      printf("%lu positive graphs: %lu vertices, %lu edges",
             parameters->numPosEgs, parameters->posGraph->numVertices,
             parameters->posGraph->numEdges);

      if (parameters->evalMethod == EVAL_MDL)
         printf(", %.0f bits\n", parameters->posGraphDL);
      else
         printf("\n");
      if (parameters->negGraph != NULL)
      {
         printf("%lu negative graphs: %lu vertices, %lu edges",
                parameters->numNegEgs, parameters->negGraph->numVertices,
                parameters->negGraph->numEdges);
         if (parameters->evalMethod == EVAL_MDL)
            printf(", %.0f bits\n", parameters->negGraphDL);
         else
            printf("\n");
      }
      printf("%lu unique labels\n", parameters->labelList->numLabels);
      printf("\n");

      // free(parameters->label_wise_count);
      // free(parameters->label_wise_mean_degree);
      parameters->label_wise_mean_degree = (double *) malloc(sizeof(double) * parameters->labelList->numLabels);
      parameters->label_wise_count = (int *) malloc(sizeof(int) * parameters->labelList->numLabels);
      for(int i=0; i<parameters->labelList->numLabels; ++i)
      {
         parameters->label_wise_count[i] = 0;
      }
      for(ULONG v=0; v<parameters->posGraph->numVertices; ++v)
      {
         parameters->posGraph->vertices[v].degree = parameters->posGraph->vertices[v].numEdges;
         parameters->label_wise_mean_degree[parameters->posGraph->vertices[v].label] += parameters->posGraph->vertices[v].degree;
         ++parameters->label_wise_count[parameters->posGraph->vertices[v].label];
      }
      for(int i=0; i<parameters->labelList->numLabels; ++i)
      {
         if(parameters->label_wise_count[i] > 0)
            parameters->label_wise_mean_degree[i] / parameters->label_wise_count[i];
      }

      subList = DiscoverSubs(parameters);
      /* Recalculate cost of each instance */
      InstanceListNode *instance_list_node = subList->head->sub->instances->head;
      double thr = parameters->threshold * (subList->head->sub->definition->numEdges + subList->head->sub->definition->numVertices);
      Graph *g_def = InstanceToGraph(instance_list_node->instance, parameters->posGraph);
      while (instance_list_node != NULL)
      {
         double match_cost;
         Graph *g = InstanceToGraph(instance_list_node->instance, parameters->posGraph);
         GraphMatch(g_def, g, parameters->labelList, thr, &match_cost, NULL, parameters->use_pose_cost);
         instance_list_node->instance->minMatchCost = match_cost;
         instance_list_node = instance_list_node->next;
      }

      /**/

      if (subList->head == NULL)
      {
         done = TRUE;
         printf("No substructures found.\n\n");
      }
      else
      {
         // printf("Calculating sub list level\n");
         // printf("Before: Label of first element: %s, level: %lu\n",
         //    parameters->labelList->labels[parameters->posGraph->vertices[0].label].labelValue.stringLabel,
         //    parameters->labelList->labels[parameters->posGraph->vertices[0].label].level);
         calculateSubListLevel(subList, parameters->labelList);
         // write output to stdout
         if (parameters->outputLevel > 3)
         {
            printf("\nBest %lu substructures:\n\n", CountSubs(subList));
            PrintSubList(subList, parameters);
            // printf("Best discovered sub cost: %f\n", subList->head->sub->instances->head->next->instance->minMatchCost);
            // EvaluateSub(subList->head->sub, parameters);
            // Graph *compressedGraph = CompressGraph(parameters->posGraph, subList->head->sub->instances, parameters);
            // SubListNode *subListNode = subList->head;
            // // instance_list_node = subList->head->sub->instances->head;
            // // while (instance_list_node != NULL)
            // // {
            // //    printf("Mean position of instance: "); printVector(instance_list_node->instance->mean_position, 3);printf("\n");
            // //    instance_list_node = instance_list_node->next;
            // // }
            // while (subListNode != NULL)
            // {
            //    printf("MDLs: compressed graph: %f, sub: %f\n", MDL(compressedGraph, parameters->labelList->numLabels, parameters), MDL(subList->head->sub->definition, parameters->labelList->numLabels, parameters));
            //    printf("Sizes: compressed graph: %f, sub: %f\n", (double)GraphSize(compressedGraph), (double)GraphSize(subListNode->sub->definition));
            //    subListNode = subListNode->next;
            // }

            /* Testing InstanceSubgraphMatch */
            // BOOLEAN *reached;
            // reached = (BOOLEAN *)malloc(sizeof(BOOLEAN) * subList->head->sub->definition->numVertices);\
            // ULONG half_pt = subList->head->sub->definition->numVertices / 2;
            // for(ULONG i=0; i<half_pt; ++i)
            // {
            //    reached[i] = FALSE;
            // }
            // for(ULONG i=half_pt; i<subList->head->sub->definition->numVertices; ++i)
            // {
            //    reached[i] = TRUE;
            // }
            // BOOLEAN match = InstanceSubgraphMatch(subList->head->sub->instances->head->instance,
            //          subList->head->sub->definition, parameters->posGraph, parameters, reached);
            // printf("Match: %d\n", match);
            /**/

            /* Testing Instance fitting */
            // markPatternInstances(subList->head->sub, parameters->posGraph, TRUE);
            // Instance *best_fitting_instance = FindFittingInstances(subList->head->sub->definition, parameters->posGraph, parameters);
            // printf("Predicted instance:\n");
            // PrintInstance(best_fitting_instance, 0, parameters->posGraph, parameters->labelList);
            // markPatternInstances(subList->head->sub, parameters->posGraph, FALSE);
            /**/
         }
         else
         {
            printf("\nBest substructure:\n\n");
            PrintSub(subList->head->sub, parameters);
         }
         printf("-1\n");

         // write machine-readable output to file, if given
         if (parameters->outputToFile)
         {
            outputFile = fopen(parameters->outFileName, "a");
            if (outputFile == NULL)
            {
               printf("WARNING: unable to write to output file %s,",
                      parameters->outFileName);
               printf("disabling\n");
               parameters->outputToFile = FALSE;
            }
            WriteGraphToFile(outputFile, subList->head->sub->definition,
                             parameters->labelList, 0, 0,
                             subList->head->sub->definition->numVertices,
                             TRUE);
            fclose(outputFile);
         }

         // Store structure and corresponding name
         // char subLabelString[TOKEN_LEN];
         // sprintf(subLabelString, "%s_%lu", SUB_LABEL_STRING, iteration);
         // strcpy(sub_names[iteration - 1], subLabelString);
         // discovered_subs[iteration - 1] = CopyGraph(subList->head->sub->definition);
         // // all_label_lists[iteration - 1] = AllocateLabelList();
         // ULONG v, e;
         // // add graph's vertex labels to new label list
         // for (v = 0; v < parameters->posGraph->numVertices; v++)
         //    parameters->posGraph->vertices[v].label =
         //        StoreLabel(&parameters->labelList->labels[parameters->posGraph->vertices[v].label],
         //                   all_label_lists[iteration - 1]);
         // // add parameters->posGraph's edge labels to new label list
         // for (e = 0; e < parameters->posGraph->numEdges; e++)
         //    parameters->posGraph->edges[e].label =
         //        StoreLabel(&parameters->labelList->labels[parameters->posGraph->edges[e].label], all_label_lists[iteration - 1]);
         // PrintGraph(discovered_subs[iteration-1], parameters->labelList);
         char subLabelString[TOKEN_LEN];
         sprintf(subLabelString, "%s_%lu", SUB_LABEL_STRING, iteration);
         // subList->head->sub->label_str = subLabelString;
         strcpy(subList->head->sub->label_str, subLabelString);
         printf("name copied \n");
         discovered_subs[iteration-1] = CopySub(subList->head->sub);
         printf("Sub copied \n");
         printf("copied: %s, struct: %s, og str: %s \n", discovered_subs[iteration-1]->label_str, subList->head->sub->label_str, subLabelString);
         printf("0\n");
         if (iteration < parameters->iterations)
         { // Another iteration?
            if (parameters->evalMethod == EVAL_SETCOVER)
            {
               printf("Removing positive examples covered by");
               printf(" best substructure.\n\n");
               RemovePosEgsCovered(subList->head->sub, parameters);
            }
            else
            {
               printf("1.1'\n");
               compressed_graphs[iteration-1] = CopyGraph(CompressFinalGraphs(subList->head->sub, parameters, iteration,
                                   FALSE));
               // for(int v=0; v<parameters->posGraph->numVertices; ++v)
               // {
               //    printf("%d (%d) : %d\n", parameters->posGraph->vertices[v].label, parameters->posGraph->vertices[v].map, compressed_graphs[iteration-1]->vertices[parameters->posGraph->vertices[v].map].label);
               //    // std::cout << parameters_->posGraph->vertices[v].label << "(" << parameters_->posGraph->vertices[v].map << ")" << " : ";
               //    // std::cout << compressed_graphs[0]->vertices[parameters_->posGraph->vertices[v].map].label << std::endl;
               // }
               compressed_graphs[parameters->iterations] = CopyGraph(parameters->posGraph);
               printf("pos graph copied\n");
               FreeGraph(parameters->posGraph);
               parameters->posGraph = CopyGraph(compressed_graphs[iteration-1]);
               all_label_lists[iteration - 1] = AllocateLabelList();
               for(ULONG li=0; li<parameters->labelList->numLabels; ++li)
               {
                  StoreLabel(&parameters->labelList->labels[li], all_label_lists[iteration-1]);
               }
            }

            // check for stopping condition
            // if set-covering, then no more positive examples
            // if MDL or size, then positive graph contains no edges
            if (parameters->evalMethod == EVAL_SETCOVER)
            {
               if (parameters->numPosEgs == 0)
               {
                  done = TRUE;
                  printf("Ending iterations - ");
                  printf("all positive examples covered.\n\n");
               }
            }
            else
            {
               if (parameters->posGraph->numEdges == 0)
               {
                  done = TRUE;
                  printf("Ending iterations - graph fully compressed.\n\n");
               }
            }
         }
         // printf("After: Label of first element: %s, level: %lu\n",
         //    parameters->labelList->labels[parameters->posGraph->vertices[0].label].labelValue.stringLabel,
         //    parameters->labelList->labels[parameters->posGraph->vertices[0].label].level);
         if ((iteration == parameters->iterations) && (parameters->compress))
         {
            if (parameters->evalMethod == EVAL_SETCOVER)
               WriteUpdatedGraphToFile(subList->head->sub, parameters);
            else
            {
               printf("2.1\n");
               compressed_graphs[iteration-1] = CopyGraph(WriteCompressedGraphToFile(subList->head->sub, parameters,
                                          iteration));
               compressed_graphs[parameters->iterations] = CopyGraph(parameters->posGraph);
               printf("pos graph copied\n");
               FreeGraph(parameters->posGraph);
               parameters->posGraph = CopyGraph(compressed_graphs[iteration-1]);
               printf("2.2\n");
               all_label_lists[iteration - 1] = AllocateLabelList();
               for(ULONG li=0; li<parameters->labelList->numLabels; ++li)
               {
                  StoreLabel(&parameters->labelList->labels[li], all_label_lists[iteration-1]);
               }
               printf("COMPRESSED GRAPH:");
               PrintGraph(compressed_graphs[iteration-1], all_label_lists[iteration-1]);
               printf("2.3\n");
            }
         }
      }
      Instance *inst = subList->head->sub->instances->head->next->instance;
      for (ULONG i = 0; i < inst->numVertices; i++)
      {
         printf("%lu -> %lu \n", inst->mapping[i].v1, inst->mapping[i].v2);
      }
      FreeSubList(subList);
      if (parameters->iterations > 1)
      {
         iterationEndTime = time(NULL);
         printf("Elapsed time for iteration %lu = %lu seconds.\n\n",
                iteration, (iterationEndTime - iterationStartTime));
      }
      iteration++;
   }

   endTime = times(&tmsend);
   printf("\nSUBDUE done (elapsed CPU time = %7.2f seconds).\n",
          (endTime - startTime) / (double)clktck);
   // for(int i=0; i<parameters->iterations; ++i)
   // {
   //    printf("Iteration %d: Sub name: %s \n", i, sub_names[i]);
   //    PrintGraph(discovered_subs[i], all_label_lists[i]);
   // }
}

Parameters *PostProcessParameters(Parameters *parameters)
{
   Increment *increment = NULL;

   if (parameters->incremental)
      increment = GetCurrentIncrement(parameters);

   // Code from this point until end of function was moved from GetParameters
   if (parameters->numPosEgs == 0)
   {
      fprintf(stderr, "ERROR: no positive graphs defined\n");
      exit(1);
   }

   // Check bounds on discovered substructures' number of vertices
   if (parameters->maxVertices == 0)
      parameters->maxVertices = parameters->posGraph->numVertices;
   if (parameters->maxVertices < parameters->minVertices)
   {
      fprintf(stderr, "ERROR: minsize exceeds maxsize\n");
      exit(1);
   }

   // Set limit accordingly
   if (parameters->limit == 0)
   {
      if (parameters->incremental)
         parameters->limit = increment->numPosEdges / 2;
      else
         parameters->limit = parameters->posGraph->numEdges / 2;
   }

   return parameters;
}