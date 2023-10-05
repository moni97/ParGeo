// This code is part of the project "ParGeo: A Library for Parallel Computational Geometry"
// Copyright (c) 2021-2022 Yiqiu Wang, Shangdi Yu, Laxman Dhulipala, Yan Gu, Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "kdTree.h"
#include "pargeo/point.h"
#include <random>
#include <cstdlib>
#include <stdlib.h>
#include <unordered_map>
#include <stdexcept>
#include <cmath>

namespace pargeo::kdTree
{
  /***********************************/
  /* Spherical range query functions */
  /***********************************/

  /********** Range Helper function **********/

  template <int dim, typename nodeT, typename objT, typename F>
  void rangeHelper(nodeT *tree, objT &q, point<dim> qMin, point<dim> qMax,
                   double radius, F func)
  {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());

    if (relation == tree->boxExclude)
    {
      return;
    }
    else if (relation == tree->boxInclude)
    {
      for (size_t i = 0; i < tree->size(); ++i)
      {
        objT *p = tree->getItem(i);
        if (p->dist(q) <= radius)
          func(p);
      }
    }
    else
    { // intersect
      if (tree->isLeaf())
      {
        for (size_t i = 0; i < tree->size(); ++i)
        {
          objT *p = tree->getItem(i);
          double dist = q.dist(*p);
          if (dist <= radius)
            func(p);
        }
      }
      else
      {
        rangeHelper<dim, nodeT, objT>(tree->L(), q, qMin, qMax, radius, func);
        rangeHelper<dim, nodeT, objT>(tree->R(), q, qMin, qMax, radius, func);
      }
    }
  }

  /********** Range Traverse function **********/

  template <int dim, typename objT, typename F>
  void rangeTraverse(
      node<dim, objT> *tree,
      objT query,
      double radius,
      F func)
  {
    point<dim> qMin, qMax;
    for (size_t i = 0; i < dim; i++)
    {
      auto tmp = query[i] - radius;
      qMin[i] = tmp;
      qMax[i] = tmp + radius * 2;
    }
    rangeHelper<dim, node<dim, objT>, objT>(tree, query, qMin, qMax,
                                            radius, func);
  }

  /********** Range Search function **********/

  template <int dim, typename objT>
  parlay::sequence<objT *> rangeSearch(
      node<dim, objT> *tree,
      objT query,
      double radius)
  {
    parlay::sequence<objT *> output;
    auto collect = [&](objT *p)
    { output.push_back(p); };
    rangeTraverse(tree, query, radius, collect);
    return output;
  }

  /********** Brute Force Range Search function **********/

  template <int dim, typename objT>
  parlay::sequence<size_t> bruteforceRange(parlay::sequence<objT> &elems,
                                           objT query,
                                           double radius)
  {
    auto out = parlay::sequence<objT>();
    auto flag = parlay::sequence<size_t>(elems.size(), elems.size());
    parallel_for(0, elems.size(), [&](size_t i)
                 {
                   if (elems[i].dist(query) <= radius)
                     flag[i] = i;
                 });
    return parlay::filter(make_slice(flag), [&](size_t i)
                          { return i < elems.size(); });
  }

  /******************************************/
  /* End of spherical range query functions */
  /******************************************/


  /************************************/
  /* Orthogonal range query functions */
  /************************************/

  /********** Orthogonal Range Helper function **********/

  template <int dim, typename nodeT, typename objT, typename F>
  void orthRangeHelper(nodeT *tree, point<dim> qMin, point<dim> qMax,
                       F func)
  {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());

    if (relation == tree->boxExclude)
    {
      return;
    }
    else if (relation == tree->boxInclude)
    {
      for (size_t i = 0; i < tree->size(); ++i)
      {
        objT *p = tree->getItem(i);
        func(p);
      }
    }
    else
    { // intersect
      if (tree->isLeaf())
      {
        for (size_t i = 0; i < tree->size(); ++i)
        {
          objT *p = tree->getItem(i);
          objT _p = *p;
          bool in = true;
          for (int d = 0; d < dim; ++d)
          {
            if (_p[d] > qMax[d] || _p[d] < qMin[d])
              in = false;
          }
          if (in)
            func(p);
        }
      }
      else
      {
        orthRangeHelper<dim, nodeT, objT>(tree->L(), qMin, qMax, func);
        orthRangeHelper<dim, nodeT, objT>(tree->R(), qMin, qMax, func);
      }
    }
  }

  /********** Orthogonal Range Traverse function **********/

  template <int dim, typename objT, typename F>
  void orthogonalRangeTraverse(node<dim, objT> *tree,
                               objT query,
                               double halfLen,
                               F func)
  {
    point<dim> qMin, qMax;
    for (size_t i = 0; i < dim; i++)
    {
      auto tmp = query[i] - halfLen;
      qMin[i] = tmp;
      qMax[i] = tmp + halfLen * 2;
    }
    orthRangeHelper<dim, node<dim, objT>, objT>(tree, qMin, qMax, func);
  }

  /********** Orthogonal Range Search function **********/

  template <int dim, typename objT>
  parlay::sequence<objT *> orthogonalRangeSearch(node<dim, objT> *tree,
                                                 objT query,
                                                 double halfLen)
  {
    parlay::sequence<objT *> output;
    auto collect = [&](objT *p)
    { output.push_back(p); };
    orthogonalRangeTraverse(tree, query, halfLen, collect);
    return output;
  }

  /********** Orthogonal Range Search Count function **********/

  template <int dim, typename objT>
  int orthogonalRangeSearchCount(node<dim, objT> *tree,
                                                 objT query,
                                                 double halfLen)
  {
    size_t pointCount = 0;
    auto collect = [&](objT *p)
    {
        ++pointCount;
    };
    orthogonalRangeTraverse(tree, query, halfLen, collect);
    return pointCount;
  }

  /********** Brute Force Orthogonal Range function **********/

  template <int dim, typename objT>
  parlay::sequence<size_t> bruteforceOrthRange(parlay::sequence<objT> &A,
                                               objT query,
                                               double halfLen)
  {
    auto out = parlay::sequence<size_t>();
    point<dim> qMin, qMax;
    for (size_t i = 0; i < dim; i++)
    {
      auto tmp = query[i] - halfLen;
      qMin[i] = tmp;
      qMax[i] = tmp + halfLen * 2;
    }
    parallel_for(0, A.size(),
                 [&](size_t i)
                 {
                   bool in = true;
                   for (int d = 0; d < dim; ++d)
                   {
                     if (A[i][d] > qMax[d] || A[i][d] < qMin[d])
                       in = false;
                   }
                   if (in)
                     out.push_back(i);
                 });
    return out;
  }

  /*******************************************/
  /* End of orthogonal range query functions */
  /******************************************/

  /*******************************************/
  /* Orthogonal range sample query functions */
  /******************************************/

  /********** Orthogonal Range Sampling Helper function **********/

  template <int dim, typename nodeT, typename objT, typename F>
  void orthSampleHelper(nodeT *tree, point<dim> qMin, point<dim> qMax,
                       F func)
  {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());
    
    if (relation == tree->boxExclude)
    {
      return;
    }
    else if (relation == tree->boxInclude)
    {
      parlay::sequence<objT*> seq(tree->items.begin(), tree->items.end());
      func(seq);
    }
    else
    { // intersect
      if (tree->isLeaf())
      {
        parlay::sequence<objT *> output;
        nodeT newNode;
        for (size_t i = 0; i < tree->size(); ++i)
        {
          objT *p = tree->getItem(i);
          objT _p = *p;
          bool in = true;
          for (int d = 0; d < dim; ++d)
          {
            if (_p[d] > qMax[d] || _p[d] < qMin[d])
              in = false;
          }
          if(in) {
            output.push_back(p);
          }
        }
        if (output.size() > 0) {
          func(output);  
        }
        // std::cout<<"Partially intersecting Node\n";
      }
      else
      {
        orthSampleHelper<dim, nodeT, objT>(tree->L(), qMin, qMax, func);
        orthSampleHelper<dim, nodeT, objT>(tree->R(), qMin, qMax, func);
      }
    }
  }

  /********** Orthogonal Range Sampling Traverse function **********/

  template <int dim, typename objT, typename F>
  void orthogonalSampleTraverse(node<dim, objT> *tree,
                               objT query,
                               double halfLen,
                               F func)
  {
    point<dim> qMin, qMax;
    for (size_t i = 0; i < dim; i++)
    {
      auto tmp = query[i] - halfLen;
      qMin[i] = tmp;
      qMax[i] = tmp + halfLen * 2;
    }
    orthSampleHelper<dim, node<dim, objT>, objT>(tree, qMin, qMax, func);
  }

  /********** Orthogonal Range Sampling function **********/

  template <int dim, typename objT>
  objT orthogonalRangeSample(node<dim, objT> *tree,
                            objT query,
                            double halfLen)
  {
    // Create a distribution (for example, uniform integer distribution)  
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Parameters
    parlay::sequence<objT*> sampleSeq;
    double M = -1.0;

    // Weighted random sampling
    // Step 1: Compute the counting query to know the number of points in Q.
    size_t numPointsInQueryRect = orthogonalRangeSearchCount(tree, query, halfLen);
    auto collect = [&](parlay::sequence<objT*> items) {
        double t = dis(gen);
        double weightOfNode = t * (double(items.size()) / numPointsInQueryRect);
        if (weightOfNode > M) {
          M = weightOfNode;
          sampleSeq = items;
        }
    };

    orthogonalSampleTraverse(tree, query, halfLen, collect);
    pargeo::point<dim> samplePoint;
    if (sampleSeq.size() > 1) {
      std::uniform_int_distribution<int> distribution(0, sampleSeq.size() - 1);
      int sampleIndex = distribution(gen);
      samplePoint = sampleSeq[sampleIndex];
    } else {
      samplePoint = sampleSeq[0];
    }
    return samplePoint;
  }

  /***************************************************/
  /* End of orthogonal range sample query functions */
  /**************************************************/

  /**************************************/
  /* Orthogonal range entropy functions */
  /**************************************/
  
  /********** Get canonical nodes function **********/

  template <int dim, typename objT>
  parlay::sequence<parlay::sequence<objT *>> getCanonicalNodes(node<dim, objT> *tree,
                            objT query,
                            double halfLen,
                            int &numCanNodes) {
    parlay::sequence<parlay::sequence<objT *>> canonicalNodes;
    auto collect = [&](parlay::sequence<objT*> items) {
        canonicalNodes.push_back(items);
        ++numCanNodes;
    };
    orthogonalSampleTraverse(tree, query, halfLen, collect);
    return canonicalNodes;
  }

  /********** Weight Sampling function **********/

  template <typename objT>
  parlay::sequence<objT *> weightedSampling(parlay::sequence<parlay::sequence<objT *>> &canonicalNodes, int numPointsInQueryRect) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    parlay::sequence<objT*> sampleSeq;
    double M = -1.0;
    double t = dis(gen);
    // Iterate through the outer sequence
    int count=0, canonicalNodesSize = 0;
    for (const auto& sequence : canonicalNodes) {
      // Iterate through the inner sequence
      double weightOfNode = t * (double(sequence.size()) / numPointsInQueryRect);
      if (weightOfNode > M) {
        M = weightOfNode;
        sampleSeq = sequence;
      }
      canonicalNodesSize += sequence.size();
      ++count;
    }
    return sampleSeq;
  }

  /********** Propotional Sampling function **********/

  template <typename objT>
  parlay::sequence<objT *> propotionalSampling(parlay::sequence<parlay::sequence<objT *>> &canonicalNodes, int numPointsInQueryRect, int m) {
    std::random_device rd;
    std::mt19937 gen(rd());

    parlay::sequence<objT*> samples;
    int i; 

    for (const auto& sequence : canonicalNodes) {
      // Iterate through the inner sequence
      double proportion = m * (double(sequence.size()) / numPointsInQueryRect);
      int numOfSamples = ceil(proportion);
      // std::cout << "seqSize: " << sequence.size() << " numOfSamples: " << numOfSamples << std::endl;
      std::uniform_int_distribution<int> distribution(0, sequence.size() - 1);
      for(i = 0; i < numOfSamples; ++i) {
        samples.push_back(sequence[distribution(gen)]);
      }  
    }
    return samples;
  }

  /********** Sample a point from a sequence function **********/
  template <typename objT>
  objT samplePointFromSeq(parlay::sequence<objT *> sequence) {
    std::random_device rd;
    std::mt19937 gen(rd());

    objT samplePoint;
    if (sequence.size() > 1) {
      std::uniform_int_distribution<int> distribution(0, sequence.size() - 1);
      int sampleIndex = distribution(gen);
      samplePoint = sequence[sampleIndex];
    } else {
      samplePoint = sequence[0];
    }
    return samplePoint;
  }

  /********** Orthogonal Range Entropy function **********/

  template <int dim, typename objT>
  double orthogonalRangeEntropyAdditive(node<dim, objT> *tree,
                    objT query,
                    double halfLen,
                    node<dim, objT> *treeMap[],
                    int n,
                    int &count,
                    int &numCanNodes,
                    double delta,
                    int m) 
  {
    /* Parameters for the algorithm */

    int i;
    count = orthogonalRangeSearchCount(tree, query, halfLen);
    
    double d, x,
      tou = (delta / n) / (10 * log2(n / delta));
    m = (log(6) / (delta * delta)) * (log2(1 / tou) * log2(1 / tou));
    
    m = fmin(count, m);
    // m = m / 10;
    double sum = 0.0, num;
    // std::cout <<  delta << ", " << m << ", ";

    /* Generates m samples, calculate d(si) and x(si) */
    // int m2 = int(m);
    std::unordered_map<int, int> numOfColorPtsInRange(n);

    parlay::sequence<parlay::sequence<objT *>> canonicalNodes = getCanonicalNodes(tree, query, halfLen, numCanNodes);

    objT s;
    for(i =0; i< n; ++i) {
      numOfColorPtsInRange[i] = orthogonalRangeSearchCount(treeMap[i], query, halfLen);
    }

    parlay::sequence<objT *> sampleSeq = propotionalSampling(canonicalNodes, count, m);
    for (i = 0; i < m; ++i) {
      s = sampleSeq[i];
      num = numOfColorPtsInRange[s.attribute];
      d =  num / count;
      /* Calculate x(i) */
      if(d >= tou)
        x = log2(double(1)/d);
      else
        x = 0;
      sum += x;
    }
    return sum / m;
  }

  template <int dim, typename objT>
  double orthogonalRangeEntropyMultiplicative(node<dim, objT> *tree,
                    objT query,
                    double halfLen,
                    node<dim, objT> *treeMap[],
                    int n,
                    int &count,
                    int &numCanNodes,
                    double delta,
                    bool isMul = false) 
  {
    /* Parameters for the algorithm */

    int i;
    count = orthogonalRangeSearchCount(tree, query, halfLen);
    
    double d, x,
      m = log2(n) / (delta * delta); 
    
    // m = fmin(count, m);
    // m = m / 100;
    double sum = 0.0, num;
    // std::cout << m << ", ";

    /* Generates m samples, calculate d(si) and x(si) */
    int m2 = int(m);
    std::unordered_map<int, int> numOfColorPtsInRange(n);
    // *numCanNodes = 0;
    parlay::sequence<parlay::sequence<objT *>> canonicalNodes = getCanonicalNodes(tree, query, halfLen, numCanNodes);

    objT s;
    for(i =0; i< n; ++i) {
      numOfColorPtsInRange[i] = orthogonalRangeSearchCount(treeMap[i], query, halfLen);
    }

    parlay::sequence<objT *> sampleSeq = propotionalSampling(canonicalNodes, count, m2);
    for (i = 0; i < m2; ++i) {
      s = sampleSeq[i];
      num = numOfColorPtsInRange[s.attribute];
      d = num / count;
      x = 0;
      if (d >= (double(1) / (n * n * n))) {
        x = log2(double(1) / d) / (3 * log2(n));
      }
      sum += x;
    }
    double res = 3 * sum * (log2(n / double(m2)));
    std::cout << delta << ", " << m << ", " << count << ", " << res << std::endl;
    return res;
  }

  /********** Orthogonal Range Entropy Brute Force function **********/

  template <int dim, typename objT>
  double rangeEntropyBruteForce(node<dim, objT> *tree,
                  objT query,
                  double halfLen,
                  int no_of_points) 
  {
    int no_of_groups = 10;
    parlay::sequence<pargeo::point<dim>*> elems2 =
        pargeo::kdTree::orthogonalRangeSearch(tree, query, halfLen);

    std::unordered_map<int, int> groupsToPointsMap(no_of_groups);

    for (pargeo::point<dim>* ptr : elems2) {
      groupsToPointsMap[ptr->attribute] += 1;
    }

    int i = 0;
    double e = 0.0;
    for(i = 0; i < no_of_groups; ++i) {
      if(groupsToPointsMap[i] != 0) {
        double div = double(groupsToPointsMap[i]) / no_of_points;
        e += div * (log2(div));
      }
    }
    return -1 * e;
  }

  /*********************************************/
  /* End of orthogonal range entropy functions */
  /*********************************************/

} // End namespace
