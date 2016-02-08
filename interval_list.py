"""Maintain a list of intervals that are merged as needed."""


class IntervalList:
    def __init__(self):
        self.ilist = []

    def __str__(self):
        return " ".join([str(x) for x in self.ilist])

    def __repr__(self):
        return str(self)

    def add(self,I):        
        """Merge a new inerval with existing overlapping intervals, or
        add it in if it doesn't overlap. 
        Optimization: if necessary, we can improve this to a binary
        search, but didn't seem worth the effort not."""

        if len(self.ilist) == 0:
            self.ilist.append(I)
            return None

        # Find: First i such that ilist[i][1] > I[0]
        a = 0
        b = len(self.ilist)

        while a < b:
            i = (a+b)//2
            if self.ilist[i][1] < I[0]:
                a = i+1
            elif i > 0 and self.ilist[i-1][1] > I[0]:
                b = i-1
            else:
                a = i
                break

        i = a
        if i == len(self.ilist):
            self.ilist.append(I)
            return None

        if I[1] <= self.ilist[i][0]:
            self.ilist.insert(i,I)
            return None

        j = i+1
        max_finish = self.ilist[i][1]
        while j < len(self.ilist) and self.ilist[j][0] < I[1]:
            max_finish = max(max_finish, self.ilist[j][1])
            j += 1

        self.ilist = self.ilist[:i] + [(min(self.ilist[i][0], I[0]), max(max_finish, I[1]))] + self.ilist[j:]
                             
    def coverage(self):
        return sum([b-a for a,b in self.ilist])

L = IntervalList()
L.add([10,20])
