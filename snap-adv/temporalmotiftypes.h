#ifndef snap_temporalmotiftypes_h
#define snap_temporalmotiftypes_h

#include "Snap.h"

typedef TKeyDat<TInt, TInt> TIntPair;

// Simple one-dimensional, two-dimensional, and three-dimensional counter
// classes with default intialization.
class Counter1D {
 public:
  Counter1D(int m=0) : m_(m) {
    if (m > 0) {
      data_ = TUInt64V(m);
      data_.PutAll(0);
    }
  }
  const TUInt64& operator()(int i) const { return data_[i]; }
  TUInt64& operator()(int i) { return data_[i]; }
  int m() { return m_; }
  
 private:
  int m_;
  TUInt64V data_;
};

class Counter2D {
 public:
  Counter2D(int m=0, int n=0) : m_(m), n_(n) {
    if (m * n > 0) {
      data_ = TUInt64V(m * n);
      data_.PutAll(0);
    }
  }
  const TUInt64& operator()(int i, int j) const { return data_[i + j * m_]; }
  TUInt64& operator()(int i, int j) { return data_[i + j * m_]; }
  int m() { return m_; }
  int n() { return n_; }
  
 private:
  int m_;
  int n_;
  TUInt64V data_;
};

class Counter3D {
 public:
  Counter3D(int m=0, int n=0, int p=0) : m_(m), n_(n), p_(p) {
    if (m * n * p > 0) {
      data_ = TUInt64V(m * n * p);
      data_.PutAll(0);
    }
  }
  const TUInt64& operator()(int i, int j, int k) const {
    return data_[i + j * n_ + k * m_ * n_];
  }
  TUInt64& operator()(int i, int j, int k) {
    return data_[i + j * m_ + k * m_ * n_];
  }
  int m() { return m_; }
  int n() { return n_; }
  int p() { return p_; }  
  
 private:
  int m_;
  int n_;
  int p_;
  TUInt64V data_;
};

// Triad edge data consists of a neighbor, a direction, and an indicator of whether
// the edge connects with wich endpoint (u or v).
class TriadEdgeData {
 public:
  TriadEdgeData() {}
  TriadEdgeData(int _nbr, int _dir, int _u_or_v) :
      nbr(_nbr), dir(_dir), u_or_v(_u_or_v) {}
  int nbr;     // Which neighbor of the center node
  int dir;     // Outgoing (0) or incoming (1) direction
  int u_or_v;  // Points to first end point u (0) or second end point v (1)
};

// Star edge data consists of a neighbor and a direction.
class StarEdgeData {
 public:
  StarEdgeData() {}
  StarEdgeData(int _nbr, int _dir) : nbr(_nbr), dir(_dir) {}
  int nbr;  // Which neighbor of the center node
  int dir;  // Outgoing (0) or incoming (1) direction
};

// TODO: THIS SHOULD BE A VIRTUAL CLASS
struct WeightFunction {
public:
  double operator()(TUInt64 start, TUInt64 end, TUInt64 window) const {
    return 0;
  }
};

struct Identity : public WeightFunction {
public:
  double operator()(TUInt64 start, TUInt64 end, TUInt64 window) const {
    return 1.;
  }
};

struct SampleWeight : public WeightFunction {
public:
  double operator()(TUInt64 start, TUInt64 end, TUInt64 window) const {
    return 1. / (1. - double(end - start) / window);
  }
};

#endif