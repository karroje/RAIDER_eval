#ifndef FAMILY_H
#define FAMILY_H
#include <vector>
#include "LmerVector.h"
#include <iostream>
#include <assert.h>

class Family {
public:
  void adopt(LmerVector *v);
  void adopt(LmerVector *v, uint L);
  
  uint size() const;
  
  uint getOffset() const { return offset; }
  void setOffset(uint off){ offset = off; }
  void addToOffset(uint off);
  
  uint getRepeatLength() const { return repeatlength; }
  void setRepeatLength(uint L){ repeatlength = L; }
  uint repeatLength(uint L) const;
  
  std::vector<LmerVector*>* getLmers() { return &vectors; }
  
  LmerVector* at(uint index) { return vectors.at(index); }
  LmerVector* see(uint index) const { return vectors.at(index); }
  void push_back(LmerVector* v) { vectors.push_back(v); }
  
  LmerVector* getPrefix() const { return vectors.front(); }
  LmerVector* getSuffix() const { return vectors.back(); }
  
  LmerVector* getLast() { return last; }
  uint getLastIndex() const { return last_index; }
  void setLast(LmerVector* v);
  bool lastRepeatComplete() const;
  
  void setExpectedEnd(uint expected) { expected_end = expected; }
  uint getExpectedEnd() const { return expected_end; }
  
  friend std::ostream &operator<<( std::ostream &output, const Family &F);
  
  LmerVector* last;
  uint last_index;
  uint expected_end;
  uint repeatlength;
  uint offset;
  std::vector<LmerVector*> vectors;
};

#endif //FAMILY_H
