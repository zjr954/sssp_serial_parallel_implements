#include <stdlib.h>       // for atoi
#include <stdio.h>        // for printf
#include <string.h>
#include <vector>
#include <queue>
#include "nodearc.h"
#include "delta.h"

using namespace std;

#define MODUL ((long long) 1 << 62)
#define VERY_FAR            9223372036854775807LL // LLONG_MAX

extern double timer(); 
extern int parse_gr(long *n_ad, long *m_ad, Node **nodes_ad, Arc **arcs_ad,
           long *node_min_ad, char *gName);
extern int parse_ss(long *sN_ad, long **source_array, char *aName);
void ArcLen(long cNodes, Node *nodes,
            long long *pMin /* = NULL */, long long *pMax /* = NULL */);
void dijkstra_serial(Node *nodes, int src_node, int n, long long *node_distance_array,
                int *visited_nodes_array, long long &cScans, long long &cUpdates);
void deltaStepping_serial(int n, int src_node, Node *nodes, int *visited_nodes_array,
		long long *tent, long long &cScans, long long &cUpdates);

int main(int argc, char **argv)
{
   Arc *arcs;
   Node *nodes; 
   long n, m, nmin, nQ;
   long long maxArcLen, minArcLen;
   long long cScans = 0LL, cUpdates = 0LL;
   long long dist;
   long *source_array=NULL;
   char gName[100], aName[100], oName[100];
   FILE *oFile;
   
   if ((argc < 4) || (argc > 5)) {
     fprintf(stderr,
             "Usage: \"%s <graph file> <aux file> <out file> [0]\"\n    or \"%s <graph file> <aux file> <out file> [<levels>] \"\n    or \"%s <graph file> <aux file> <out file> [-<log delta>] \"\n", argv[0], argv[0], argv[0]);
     exit(0);
   }
   strcpy(gName, argv[1]);
   strcpy(aName, argv[2]);
   strcpy(oName, argv[3]);
   oFile = fopen(oName, "a");
  
   fprintf(stderr,"c ---------------------------------------------------\n");
   fprintf(stderr,"c SSSP serial version by Xutong Ren\n");
   fprintf(stderr,"c ---------------------------------------------------\n");

   parse_gr(&n, &m, &nodes, &arcs, &nmin, gName);
   parse_ss(&nQ, &source_array, aName);
   ArcLen(n, nodes, &minArcLen, &maxArcLen);

   fprintf(oFile, "f %s %s\n", gName, aName);
   fprintf(stderr,"c\n");
   fprintf(stderr,"c Nodes: %24ld       Arcs: %22ld\n",  n, m);
   fprintf(stderr,"c MinArcLen: %20lld       MaxArcLen: %17lld\n",
           minArcLen, maxArcLen);
   fprintf(stderr,"c Trials: %23ld\n", nQ);

   long long * node_distance_array = (long long *) malloc (n * sizeof(long long));
   int * visited_nodes_array = (int *) malloc (n * sizeof(int));
   double tm = 0.0;
   tm = timer();          // start timing

   for (int i = 0; i < nQ; i++) {
       dijkstra_serial(nodes, source_array[i] - 1, n, node_distance_array,
		visited_nodes_array, cScans, cUpdates);
       //deltaStepping_serial(n, source_array[i] - 1, nodes, visited_nodes_array,
	//		node_distance_array, cScans, cUpdates);
#ifdef CHECKSUM
       dist = 0LL;
       for (long i = 0; i < n; ++i)
         if (visited_nodes_array[i]) {
           dist = (dist + (node_distance_array[i] % MODUL)) % MODUL;
         }
       fprintf(oFile,"d %lld\n", dist);
#endif
   }
   tm = (timer() - tm);   // finish timing

#ifndef CHECKSUM
   fprintf(stderr, "c Scans (ave): %20.1f     Improvements (ave): %10.1f\n",
          (float) cScans / (float) nQ,
          (float) cUpdates / (float) nQ);
   fprintf(stderr,"c Time (ave, ms): %18.2f\n",
            1000.0 * tm/(float) nQ);
   fprintf(oFile, "g %ld %ld %lld %lld\n",
           n, m, minArcLen, maxArcLen);
   fprintf(oFile, "t %f\n", 1000.0 * tm/(float) nQ);
   fprintf(oFile, "v %f\n", (float) cScans/ (float) nQ);
   fprintf(oFile, "i %f\n", (float) cUpdates/ (float) nQ);
#endif

   free(node_distance_array);
   free(visited_nodes_array);
   free(source_array);
   fclose(oFile);
   return 0;
}

void ArcLen(long cNodes, Node *nodes,
            long long *pMin /* = NULL */, long long *pMax /* = NULL */)
{
   Arc *lastArc, *arc;
   long long maxLen = 0, minLen = VERY_FAR;
   lastArc = (nodes+cNodes)->first - 1;
   for ( arc = nodes->first; arc <= lastArc; arc++ )
   {
      if ( arc->len > maxLen )
         maxLen = arc->len;
      if ( arc->len < minLen )
         minLen = arc->len;
   }
   if ( pMin )   *pMin = minLen;
   if ( pMax )   *pMax = maxLen;
}

/*************************** Serial implementation for Dijkstra ***************************/
/* Dijkstra with priority queue */
void dijkstra_serial(Node *nodes, int src_node, int n, long long *node_distance_array,
                int *visited_nodes_array, long long &cScans, long long &cUpdates)
{
  /* Initialization */
  for (long i = 0; i < n; i++) {
    node_distance_array[i] = VERY_FAR;
    visited_nodes_array[i] = 0;
  }
  node_distance_array[src_node] = 0;
  priority_queue<Arc> q;
  Arc s;
  s.len = 0LL, s.head = nodes + src_node;
  q.push(s);

  /* Dijkstra Core */
  while(!q.empty()) {
    s = q.top(), q.pop();
    if (visited_nodes_array[s.head - nodes]) continue;
    long cur_node = s.head - nodes;
    visited_nodes_array[cur_node] = 1;
    cScans++;
    Arc *lastArc, *arc;
    arc = s.head->first;
    lastArc = (s.head + 1)->first;
    for (; arc < lastArc; arc++) {
      long next = arc->head - nodes;
      if (visited_nodes_array[next] != 1) {
        long long new_distance = node_distance_array[cur_node] + arc->len;
        /* Relaxation */
        if (new_distance < node_distance_array[next]) {
          node_distance_array[next] = new_distance;
          s.len = new_distance, s.head = arc->head;
          q.push(s);
          cUpdates++;
        }
      } //end if
    } //end for arcs
  } //end for nodes
}
/*************************** Serial implementation for Dijkstra ***************************/

/*************************** Serial implementation for Delta-Stepping ***************************/
void relax(int *buckets, int *bCnt, long long *tent, struct request req, long long &cUpdates)
{
        if (req.weight < tent[req.vertex]) {
                cUpdates++;
                if (tent[req.vertex] != VERY_FAR)
                {
                        bCnt[buckets[req.vertex]]--;
                        if (bCnt[buckets[req.vertex]] >= 0);
                        buckets[req.vertex] = -1;
                }
                int index = ((int)((double)req.weight / delta)) % bucketsCount;
                if (buckets[req.vertex] != index) {
                        if (buckets[req.vertex] != -1)
                                bCnt[buckets[req.vertex]]--;
                        buckets[req.vertex] = index;
                        bCnt[index]++;
                }
                tent[req.vertex] = req.weight;
        }
}

void findRequests(int n, Node* nodes, long long *tent, int *bucket, int minIndex, KIND kind, vector<struct request> &requests, long long &cScans)
{
        for (int it = 0; it < n; ++it)
        {
                if (bucket[it] != minIndex) continue;
                cScans++;
                Arc *arc = (nodes+ it)->first;
                Arc *lastArc = (nodes+ it + 1)->first;
                if (kind == LIGHT)
                {
                        for (Arc *i = arc; i < lastArc; i++)
                        {
                                if (i->len <= delta)
                                {
                                        struct request req;
                                        req.vertex = i->head - nodes;
                                        req.weight = tent[it] + i->len;
                                        requests.push_back(req);
                                }
                        }
                }
                else
                {
                        for (Arc *i = lastArc - 1; i >= arc; i--)
                        {
                                if (i->len > delta)
                                {
                                        struct request req;
                                        req.vertex = i->head - nodes;
                                        req.weight = tent[it] + i->len;
                                        requests.push_back(req);
                                }
                        }
                }
        }
}

void deltaStepping_serial(int n, int source, Node *nodes, int * buckets, long long *tent, long long &cScans, long long &cUpdates)
{
        /* Initialization */
        for (int i = 0; i < n; i++) {
            tent[i] = VERY_FAR;
        }
        int bCnt[bucketsCount];
        memset(buckets, -1, n*sizeof(int));
        memset(bCnt, 0, bucketsCount*sizeof(int));
        buckets[source] = 0;
        bCnt[0] = 1;
        tent[source] = 0;
        int cyclingIndex = 0;

        while (true)
        {
                int minIndex = -1;
                for (int i = 0; i < bucketsCount; i++)
                {
                        int indexForCheck = i + cyclingIndex;
                        if (indexForCheck >= bucketsCount)
                        {
                                indexForCheck -= bucketsCount;
                        }
                        if (bCnt[indexForCheck] > 0)
                        {
                                minIndex = i;
                                break;
                        }
                }
                cyclingIndex++;
                if (cyclingIndex >= bucketsCount)
                        cyclingIndex = 0;

                if (minIndex == -1)
                {
                        break;
                }

                int currentSet[n];
                memset(currentSet, -1, n*sizeof(int));

                while(bCnt[minIndex])
                {
                        vector<struct request> requests;
                        findRequests(n, nodes, tent, buckets, minIndex, LIGHT, requests, cScans);

                        for (int i = 0; i < n; ++i)
                        {
                                if (buckets[i] == minIndex) {
                                        currentSet[i] = minIndex;
                                        buckets[i] = -1;
                                }
                        }
                        bCnt[minIndex] = 0;

                        for (int i = 0; i < requests.size(); i++)
                        {
                                relax(buckets, bCnt, tent, requests[i], cUpdates);
                        }
                        requests.clear();
                }

                vector<struct request> heavyRequests;
                findRequests(n, nodes, tent, currentSet, minIndex, HEAVY, heavyRequests, cScans);
                for (int i = 0; i < heavyRequests.size(); i++)
                {
                        relax(buckets, bCnt, tent, heavyRequests[i], cUpdates);
                }
                heavyRequests.clear();
        } // end while
}
/*************************** Serial implementation for Delta-Stepping ***************************/
