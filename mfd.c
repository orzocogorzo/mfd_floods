#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>



int * drainpaths (
    float start,
    float flow,
    int xsize,
    int ysize
) {
    static int next_step[10];
    
    for rc 
}

void main (self, start, flow):
    self.drainages = np.zeros(self.dtm.shape)
    
    def _drainpaths (rcs, step_drainages, queue=list(), visited=dict()):
        try:
            next_step = list()

            for rc in rcs:
                if rc in visited: continue                        
                visited[rc] = True
                    
                src_flood = step_drainages[rc]
                if src_flood < 1e-1:
                    src_flood += self.drainages[rc]
                    if src_flood < 1e-2: continue
                
                perimetters = self.get_perimetters(*rc)
                downslopes = self.get_downslopes(perimetters)
                if len(downslopes) == 0:
                    # LOCAL DEPRESSION SCENARIO
                    src_flood = self.drainages[rc] - self.drainages[(
                        self.deltas[perimetters.argmin()][0],
                        self.deltas[perimetters.argmin()][1]
                    )] + src_flood
                    overflow = src_flood - self.get_volumetry(perimetters)
                    if overflow > 0:
                        delta = self.deltas[perimetters.argmin()]
                        new_rc = (rc[0] + delta[0], rc[1] + delta[1])
                        step_drainages[new_rc] += overflow
                        if new_rc not in visited: next_step.append(new_rc)
                else:
                    t_downslopes = reduce(lambda a, d: a + min(0, d), downslopes.values())
                    if t_downslopes == 0:
                        # ABSOLUTE PLANE SCENARIO
                        for delta in downslopes:
                            new_rc = (rc[0] + delta[0], rc[1] + delta[1])
                            catchment_factor = 1/len(downslopes)
                            new_catchment = catchment_factor * src_flood
                            step_drainages[new_rc] += new_catchment
                            if new_rc not in visited: next_step.append(new_rc)
                    else:
                        # DOWNSLOPES SCENARIO
                        for delta in downslopes:
                            new_rc = (rc[0] + delta[0], rc[1] + delta[1])
                            catchment_factor = downslopes[delta] / t_downslopes
                            new_catchment = catchment_factor * src_flood
                            step_drainages[new_rc] += new_catchment
                            if new_rc not in visited: next_step.append(new_rc)
            
            if len(next_step): queue.append(next_step)
            if len(queue) > 0: _drainpaths(queue.pop(), step_drainages, queue, visited)
            
            return step_drainages
                    
        except RuntimeWarning as e:
            print("RuntimeWarning")
            print(e)
        except KeyboardInterrupt as e:
            print("KeyboardInterrupt!!")
            return step_drainages
        except ValueError as e:
            print("ValueError!!")
            print(e)
            return step_drainages
        except RecursionError as e:
            print("RecursionError!!")
            print(e)
            return step_drainages
        except Exception as e:
            print("Exception!!")
            print(e)
            return step_drainages

    # print(flow)
    # self.drainages[start] = float(flow)
    # self.drainages = _drainpaths([start], self.drainages)  
    hyd = hydrogram(flow)
    for flood in hyd:
        step_drainages = np.zeros(self.drainages.shape)
        step_drainages[start] = float(flood)
        step_drainages = _drainpaths([start], step_drainages, list(), dict())
        self.drainages = np.where(self.drainages >= step_drainages, self.drainages, step_drainages)

    return self.drainages


