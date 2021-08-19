/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2015-2017
 * EMAXA Kutato-fejleszto Kft. (EMAXA Research Ltd.)
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#include <lemon/vf2.h>
#include <lemon/vf2pp.h>
#include <lemon/concepts/digraph.h>
#include <lemon/smart_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/concepts/maps.h>


#include "test/test_tools.h"
#include <sstream>

using namespace lemon;

char petersen_lgf[] =
  "@nodes\n"
  "label col1 col2\n"
  "0 1 1\n"
  "1 1 2\n"
  "2 1 3\n"
  "3 1 4\n"
  "4 2 5\n"
  "5 2 1\n"
  "6 2 2\n"
  "7 2 3\n"
  "8 2 4\n"
  "9 2 5\n"
  "@arcs\n"
  "     -\n"
  "0 1\n"
  "1 2\n"
  "2 3\n"
  "3 4\n"
  "4 0\n"
  "0 5\n"
  "1 6\n"
  "2 7\n"
  "3 8\n"
  "4 9\n"
  "5 8\n"
  "5 7\n"
  "9 6\n"
  "9 7\n"
  "6 8\n";

char c5_lgf[] =
  "@nodes\n"
  "label col\n"
  "0 1\n"
  "1 2\n"
  "2 3\n"
  "3 4\n"
  "4 5\n"
  "@arcs\n"
  "     -\n"
  "0 1\n"
  "1 2\n"
  "2 3\n"
  "3 4\n"
  "4 0\n";

char c7_lgf[] =
  "@nodes\n"
  "label\n"
  "0\n"
  "1\n"
  "2\n"
  "3\n"
  "4\n"
  "5\n"
  "6\n"
  "@arcs\n"
  "     -\n"
  "0 1\n"
  "1 2\n"
  "2 3\n"
  "3 4\n"
  "4 5\n"
  "5 6\n"
  "6 0\n";

char c10_lgf[] =
  "@nodes\n"
  "label\n"
  "0\n"
  "1\n"
  "2\n"
  "3\n"
  "4\n"
  "5\n"
  "6\n"
  "7\n"
  "8\n"
  "9\n"
  "@arcs\n"
  "     -\n"
  "0 1\n"
  "1 2\n"
  "2 3\n"
  "3 4\n"
  "4 5\n"
  "5 6\n"
  "6 7\n"
  "7 8\n"
  "8 9\n"
  "9 0\n";

char p10_lgf[] =
  "@nodes\n"
  "label\n"
  "0\n"
  "1\n"
  "2\n"
  "3\n"
  "4\n"
  "5\n"
  "6\n"
  "7\n"
  "8\n"
  "9\n"
  "@arcs\n"
  "     -\n"
  "0 1\n"
  "1 2\n"
  "2 3\n"
  "3 4\n"
  "4 5\n"
  "5 6\n"
  "6 7\n"
  "7 8\n"
  "8 9\n";

SmartGraph petersen, c5, c7, c10, p10;
SmartGraph::NodeMap<int> petersen_col1(petersen);
SmartGraph::NodeMap<int> petersen_col2(petersen);
SmartGraph::NodeMap<int> c5_col(c5);

void make_graphs() {
  std::stringstream ss(petersen_lgf);
  graphReader(petersen, ss)
    .nodeMap("col1", petersen_col1)
    .nodeMap("col2", petersen_col2)
    .run();

  ss.clear();
  ss.str("");
  ss << c5_lgf;
  //std::stringstream ss2(c5_lgf);
  graphReader(c5, ss)
    .nodeMap("col", c5_col)
    .run();

  ss.clear();
  ss.str("");
  ss << c7_lgf;
  graphReader(c7, ss).run();

  ss.clear();
  ss.str("");
  ss << c10_lgf;
  graphReader(c10, ss).run();

  ss.clear();
  ss.str("");
  ss << p10_lgf;
  graphReader(p10, ss).run();

}

class EqComparable {
public:
  bool operator==(const EqComparable&) {
    return false;
  }
};

template<class A, class B>
class EqClass {
public:
  bool operator()(A, B){
    return false;
  }
};

class IntConvertible1 {
public:
  operator int() {
    return 0;
  }
};

class IntConvertible2 {
public:
  operator int() {
    return 0;
  }
};

template<class G1,class G2>
void checkVf2Compile() {
  G1 g;
  G2 h;
  concepts::ReadWriteMap<typename G1::Node, typename G2::Node> r;
  bool succ;
  ::lemon::ignore_unused_variable_warning(succ);

  succ = vf2(g,h).run();
  succ = vf2(g,h).induced().run();
  succ = vf2(g,h).iso().run();
  succ = vf2(g,h).mapping(r).run();

  Vf2<G1,G2,concepts::ReadWriteMap<typename G1::Node, typename G2::Node>,
      EqClass<typename G1::Node,typename G2::Node> >
    myVf2(g,h,r,EqClass<typename G1::Node,typename G2::Node>());
  myVf2.find();

  succ = vf2(g,h).induced().mapping(r).run();
  succ = vf2(g,h).iso().mapping(r).run();

  concepts::ReadMap<typename G1::Node, EqComparable> l1;
  concepts::ReadMap<typename G2::Node, EqComparable> l2;
  succ = vf2(g,h).nodeLabels(l1,l2).mapping(r).run();
  succ = vf2(g,h).nodeEq(EqClass<typename G1::Node,typename G2::Node>())
    .mapping(r).run();
}

void vf2Compile() {
  checkVf2Compile<concepts::Graph,concepts::Graph>();
  checkVf2Compile<concepts::Graph,SmartGraph>();
  checkVf2Compile<SmartGraph,concepts::Graph>();
}

template<class G1,class G2>
void checkVf2ppCompile() {
  G1 g;
  G2 h;
  concepts::ReadWriteMap<typename G1::Node, typename G2::Node> r;
  bool succ;
  ::lemon::ignore_unused_variable_warning(succ);

  succ = vf2pp(g,h).run();
  succ = vf2pp(g,h).induced().run();
  succ = vf2pp(g,h).iso().run();
  succ = vf2pp(g,h).mapping(r).run();
  succ = vf2pp(g,h).induced().mapping(r).run();
  succ = vf2pp(g,h).iso().mapping(r).run();

  concepts::ReadMap<typename G1::Node, int> c1;
  concepts::ReadMap<typename G2::Node, int> c2;
  Vf2pp<G1,G2,concepts::ReadWriteMap<typename G1::Node, typename G2::Node>,
        concepts::ReadMap<typename G1::Node, int>,
        concepts::ReadMap<typename G2::Node, int> >
    myVf2pp(g,h,r,c1,c2);
  myVf2pp.find();

  succ = vf2pp(g,h).nodeLabels(c1,c2).mapping(r).run();
  succ = vf2pp(g,h).nodeLabels(c1,c2).mapping(r).run();

  concepts::ReadMap<typename G1::Node, char> c1_c;
  concepts::ReadMap<typename G2::Node, char> c2_c;
  Vf2pp<G1,G2,concepts::ReadWriteMap<typename G1::Node, typename G2::Node>,
        concepts::ReadMap<typename G1::Node, char>,
        concepts::ReadMap<typename G2::Node, char> >
    myVf2pp_c(g,h,r,c1_c,c2_c);
  myVf2pp_c.find();

  succ = vf2pp(g,h).nodeLabels(c1_c,c2_c).mapping(r).run();
  succ = vf2pp(g,h).nodeLabels(c1_c,c2_c).mapping(r).run();

  concepts::ReadMap<typename G1::Node, IntConvertible1> c1_IntConv;
  concepts::ReadMap<typename G2::Node, IntConvertible2> c2_IntConv;
  Vf2pp<G1,G2,concepts::ReadWriteMap<typename G1::Node, typename G2::Node>,
        concepts::ReadMap<typename G1::Node, IntConvertible1>,
        concepts::ReadMap<typename G2::Node, IntConvertible2> >
    myVf2pp_IntConv(g,h,r,c1_IntConv,c2_IntConv);
  myVf2pp_IntConv.find();

  succ = vf2pp(g,h).nodeLabels(c1_IntConv,c2_IntConv).mapping(r).run();
  succ = vf2pp(g,h).nodeLabels(c1_IntConv,c2_IntConv).mapping(r).run();
}

void vf2ppCompile() {
  checkVf2ppCompile<concepts::Graph,concepts::Graph>();
  checkVf2ppCompile<concepts::Graph,SmartGraph>();
  checkVf2ppCompile<SmartGraph,concepts::Graph>();
}

template<class G1, class G2, class I>
void checkSubIso(const G1 &g1, const G2 &g2, const I &i) {
  std::set<typename G2::Node> image;
  for (typename G1::NodeIt n(g1);n!=INVALID;++n){
    check(i[n]!=INVALID, "Wrong isomorphism: incomplete mapping.");
    check(image.count(i[n])==0,"Wrong isomorphism: not injective.");
    image.insert(i[n]);
  }
  for (typename G1::EdgeIt e(g1);e!=INVALID;++e) {
    check(findEdge(g2,i[g1.u(e)],i[g1.v(e)])!=INVALID,
          "Wrong isomorphism: missing edge(checkSub).");
  }
}

template<class G1, class G2, class I>
void checkIndIso(const G1 &g1, const G2 &g2, const I &i) {
  std::set<typename G2::Node> image;
  for (typename G1::NodeIt n(g1);n!=INVALID;++n) {
    check(i[n]!=INVALID, "Wrong isomorphism: incomplete mapping.");
    check(image.count(i[n])==0,"Wrong isomorphism: not injective.");
    image.insert(i[n]);
  }
  for (typename G1::NodeIt n(g1); n!=INVALID; ++n) {
    for (typename G1::NodeIt m(g1); m!=INVALID; ++m) {
      if((findEdge(g1,n,m)==INVALID) != (findEdge(g2,i[n],i[m])==INVALID)) {
        std::cout << "Wrong isomorphism: edge mismatch";
        exit(1);
      }
    }
  }
}

template<class G1,class G2,class T>
bool checkSub(const G1 &g1, const G2 &g2, const T &vf2) {
  typename G1:: template NodeMap<typename G2::Node> iso(g1,INVALID);
  if (const_cast<T&>(vf2).mapping(iso).run()) {
    checkSubIso(g1,g2,iso);
    return true;
  }
  return false;
}

template<class G1,class G2,class T>
bool checkInd(const G1 &g1, const G2 &g2, const T &vf2) {
  typename G1:: template NodeMap<typename G2::Node> iso(g1,INVALID);
  if (const_cast<T&>(vf2).induced().mapping(iso).run()) {
    checkIndIso(g1,g2,iso);
    return true;
  }
  return false;
}

template<class G1,class G2,class T>
bool checkIso(const G1 &g1, const G2 &g2, const T &vf2) {
  typename G1:: template NodeMap<typename G2::Node> iso(g1,INVALID);
  if (const_cast<T&>(vf2).iso().mapping(iso).run()) {
    check(countNodes(g1)==countNodes(g2),
          "Wrong iso alg.: they are not isomophic.");
    checkIndIso(g1,g2,iso);
    return true;
  }
  return false;
}

template<class G1, class G2, class L1, class L2, class I>
void checkLabel(const G1 &g1, const G2 &,
                const L1 &l1, const L2 &l2, const I &i) {
  for (typename G1::NodeIt n(g1);n!=INVALID;++n) {
    check(l1[n]==l2[i[n]],"Wrong isomorphism: label mismatch.");
  }
}

template<class G1,class G2,class L1,class L2,class T>
bool checkSub(const G1 &g1, const G2 &g2, const L1 &l1, const L2 &l2, const T &vf2) {
  typename G1:: template NodeMap<typename G2::Node> iso(g1,INVALID);
  if (const_cast<T&>(vf2).nodeLabels(l1,l2).mapping(iso).run()){
    checkSubIso(g1,g2,iso);
    checkLabel(g1,g2,l1,l2,iso);
    return true;
  }
  return false;
}

template<class G1,class G2>
void checkSub(const G1 &g1, const G2 &g2, bool expected, const char* msg) {
  check(checkSub(g1,g2,vf2(g1,g2)) == expected, msg);
  check(checkSub(g1,g2,vf2pp(g1,g2)) == expected, msg);
}

template<class G1,class G2>
void checkInd(const G1 &g1, const G2 &g2, bool expected, const char* msg) {
  check(checkInd(g1,g2,vf2(g1,g2)) == expected, msg);
  check(checkInd(g1,g2,vf2pp(g1,g2)) == expected, msg);
}

template<class G1,class G2>
void checkIso(const G1 &g1, const G2 &g2, bool expected, const char* msg) {
  check(checkIso(g1,g2,vf2(g1,g2)) == expected, msg);
  check(checkIso(g1,g2,vf2pp(g1,g2)) == expected, msg);
}

template<class G1,class G2,class L1,class L2>
void checkSub(const G1 &g1, const G2 &g2, const L1 &l1, const L2 &l2, bool expected, const char* msg) {
  check(checkSub(g1,g2,l1,l2,vf2(g1,g2)) == expected, msg);
  check(checkSub(g1,g2,l1,l2,vf2pp(g1,g2)) == expected, msg);
}

int main() {
	lemon::ListGraph pattern = ListGraph();
	auto n1 =pattern.addNode();
	auto n2 =pattern.addNode();
	auto n3 =pattern.addNode();
	pattern.addEdge(n1,n2);
	pattern.addEdge(n2,n3);

	ListGraph target = ListGraph();
	auto nn1 =target.addNode();
	auto nn2 =target.addNode();
	auto nn3 =target.addNode();
	auto nn4 =target.addNode();
	auto nn5 =target.addNode();
	target.addEdge(nn1,nn2);
	target.addEdge(nn2,nn3);
	target.addEdge(nn1,nn3);
	target.addEdge(nn3,nn4);
	target.addEdge(nn3,nn5);
	target.addEdge(nn4,nn5);

	int num_isos = vf2pp(pattern,target).iso().count();
	std::cout << num_isos << '\n';

    ListGraph::NodeMap<ListGraph::Node> m(pattern);
	ListGraph::NodeMap<int> map1(pattern);
	map1[n1] = 1;
	map1[n2] = 1;
	map1[n3] = 2;
	ListGraph::NodeMap<int> map2(target);
	map2[nn1] =1;
	map2[nn2] =1;
	map2[nn3] =2;
	map2[nn4] =2;
	map2[nn5] =2;
	ListGraph::NodeMap<int> map3(target);
	map3[nn1] =1;
	map3[nn2] =2;
	map3[nn3] =3;
	map3[nn4] =4;
	map3[nn5] =5;
	auto* myVf2pp = vf2pp(pattern,target).mapping(m).nodeLabels(map1,map2).getPtrToVf2ppObject();
	int c = 0;
    while(myVf2pp->find()){
      //process the current mapping m
	  c++;
	  std::cout << map3[m[n1]] << " " << map3[m[n2]] << " " << map3[m[n3]] << '\n';
    }
	delete myVf2pp;
	std::cout << c << '\n';
}
