#pragma once
#include <utility>
#include <memory>
#include <functional>
#include <algorithm>

#include "conan/mvector.h"

/*
 * KdTree algorithm
 * To count the number of points that satisfy a certain
 * predicate, you have to define a function object
 * that derives from predicate.
 * This object needs not be constant. So that a comparison
 * with a bounding box can influence the way points within 
 * this box are tested.
 */

namespace kdTree {

// should change to some different template scheme

template <typename Point, int R>
class BoundingBox
{
	std::function<double (Point const &, int)> index;

	public:
		Conan::mVector<double, R> X1, X2;
		BoundingBox(std::function<double (Point const &, int)> index_):
       			index(index_) {}

		BoundingBox(
			Conan::mVector<double, R> const &X1_, 
			Conan::mVector<double, R> const &X2_,
			std::function<double (Point const &, int)> index_);

		template <typename Iter>
		BoundingBox(Iter begin, Iter end,
			std::function<double (Point const &, int)> index_);

		double max_coord(unsigned i) const { return X2[i]; }
		double min_coord(unsigned i) const { return X1[i]; }
		double half(unsigned i);
};

template <typename Point, int R>
class Distance
{
	public:
		virtual double operator()(Point const &) const = 0;
		virtual double operator()(BoundingBox<Point, R> const &) const = 0;
		virtual ~Distance() {}
};

template <typename Point, int R>
class Predicate
{
	public:
		virtual bool operator()(Point const &) const = 0;
		virtual bool operator()(BoundingBox<Point, R> const &) const = 0;
		virtual ~Predicate() {}
};

template <typename Point>
class Visitor
{
	public:
		virtual void operator()(Point const &) = 0;
		virtual ~Visitor() {}
};

template <typename Point, int R>
class Branch: public BoundingBox<Point, R>
{
	public:
		Branch() {}

		template <typename Iter>
		Branch(Iter begin, Iter end,
			std::function<double (Point const &, int)> idx): 
			BoundingBox<Point, R>(begin, end, idx) {}

		virtual size_t count_if(Predicate<Point, R> const &pred) const = 0;
		virtual ~Branch() {}
		virtual void traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred) const = 0;
		virtual Point const &nearest_neighbour(Distance<Point, R> const &dist) const = 0;
};

template <typename Point, int R>
class Tree: public Branch<Point, R>
{
	typedef std::unique_ptr<Branch<Point, R>> Branch_ptr;
	std::pair<Branch_ptr, Branch_ptr> B;

	public:
		template <typename Iter>
		Tree(Iter begin, Iter end, 
			std::function<double (Point const &, int)> idx, 
			int dim = 0, int depth = 0);

		virtual size_t count_if(Predicate<Point, R> const &pred) const;
		virtual void traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred) const;
		virtual Point const &nearest_neighbour(Distance<Point, R> const &dist) const;
//	private:
//		Tree(Tree const &); //NI
};

template <typename Point, int R, typename Iter>
class Leaf: public Branch<Point, R>
{
	std::pair<Iter, Iter>	I;

	public:
		Leaf(Iter begin, Iter end,
			std::function<double (Point const &, int)> idx);

		virtual size_t count_if(Predicate<Point, R> const &pred) const;
		virtual void traverse(Visitor<Point> &visit, Predicate<Point, R> const &pred) const;
		virtual Point const &nearest_neighbour(Distance<Point, R> const &dist) const;
};

} // namespace kdTree

#include "kdtree_implementation.h"

