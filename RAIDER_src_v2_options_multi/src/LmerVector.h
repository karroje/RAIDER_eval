#ifndef LMER_VECTOR_H
#define LMER_VECTOR_H
#include <vector>
#include <stdexcept>
#include <iostream>

#ifndef uint
typedef unsigned int uint;
#endif

class Family;

class LmerVector {
public:
  
  LmerVector(uint position);
  ~LmerVector() { }
  
  void push_back(uint val);
  
  uint operator [](uint i) const { return lmers[i]; }
  uint & operator [](uint i) { return lmers[i]; }
  
  uint prev() const { return previous; }
  uint back() const { return lmers.back(); }
  uint front() const { return lmers.front(); }
  
  uint getPosition() const { return position; }
  void setPosition(uint i) { position = i; }
  
  uint size() const { return lmers.size(); }
  
  Family* getFamily() const { return family; }
  void setFamily(Family *fam) { family = fam; }
  
  friend std::ostream &operator<<( std::ostream &output, const LmerVector &L);
  
  std::vector<uint> lmers;
  Family* family;
  uint previous;
  uint position;
};

#endif //LMER_VECTOR_H
