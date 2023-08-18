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

namespace pargeo::kdTree
{

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

  /* Sampling in a range */

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
    node<dim, objT>* sampleNodePtr = nullptr;
    double M = -1.0;

    // Weighted random sampling

    // Step 1: Compute the counting query to know the number of points in Q.
    size_t numPointsInQueryRect = orthogonalRangeSearchCount(tree, query, halfLen);
    // std::cout << "num of points: " << numPointsInQueryRect << std::endl;
    auto collect = [&](node<dim, objT>* n) {
        double t = dis(gen);
        double weightOfNode = t * (double(n->size()) / numPointsInQueryRect);
        // std::cout << "n size: " << n->size() << " t : " << t << " weightOfNode: " << weightOfNode << " n / w: " <<double(n->size()) / numPointsInQueryRect<< std::endl;
        if (weightOfNode > M) {
          M = weightOfNode;
          sampleNodePtr = n;
        }
    };

    orthogonalSampleTraverse(tree, query, halfLen, collect);
    pargeo::point<dim> samplePoint;
    if (sampleNodePtr->size() > 1) {
      std::uniform_int_distribution<int> distribution(1, sampleNodePtr->size());
      int sampleIndex = distribution(gen);
      samplePoint = sampleNodePtr->getItem(sampleIndex);
    } else {
      samplePoint = sampleNodePtr->getItem(0);
    }

    return samplePoint;
  }

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
      
      func(tree);
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
        newNode = nodeT();
        newNode.items = parlay::make_slice(output);
        func(&newNode);
      }
      else
      {
        orthSampleHelper<dim, nodeT, objT>(tree->L(), qMin, qMax, func);
        orthSampleHelper<dim, nodeT, objT>(tree->R(), qMin, qMax, func);
      }
    }
  }

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

  template <int dim, typename objT>
  double orthogonalRangeEntropy(node<dim, objT> *tree,
                    objT query,
                    double halfLen,
                    node<dim, objT> *treeMap[]) 
    {
        /* Parameters for the algorithm */

        int i, 
          n = 3, 
          count = orthogonalRangeSearchCount(tree, query, halfLen);
      
        double d,
          delta = 0.1,
          tou = (delta / n) / (10 * log2(n / delta)),
          m = (log(6) / (delta * delta)) * (log2(1 / tou) * log2(1 / tou));
          
        /* since m is too large we take the minimum of m and count */
        
        m = fmin(count, m);

        double x[int(m)], 
          sum = 0.0;
        point<dim> s[int(m)];

        /* print all the values */

        std::cout << "\nn : " << n << "\ndelta : " << delta << "\ntou : " << tou << "\nm : " << m << std::endl;
        
        /* Generates m samples, calculate d(si) and x(si) */

        for (i = 0; i < m; ++i) {
          s[i] = orthogonalRangeSample(tree, query, halfLen);
          d =  double(orthogonalRangeSearchCount(treeMap[s[i].attribute], query, halfLen)) / count;
    
          /* Calculate x(i) */
          if(d >= tou)
            x[i] = log(double(1)/d);
          else
            x[i] = 0;
          sum += x[i];
        }

        return sum / m;
    }

} // End namespace
