def bboxtraj(atraj, margin = 0.):
    bbox = []
    bbox.append(min(atraj.data['geometry'].x.values) - margin)
    bbox.append(max(atraj.data['geometry'].x.values) + margin)
    bbox.append(min(atraj.data['geometry'].y.values) - margin)
    bbox.append(max(atraj.data['geometry'].y.values) + margin)
    return(bbox)

def bboxtrajgroup(agroup, margin = 0.):
    bbox = []
    allx = np.concatenate([ p.data['geometry'].x.values[:] for p in agroup ])
    ally = np.concatenate([ p.data['geometry'].y.values[:] for p in agroup ])
    bbox.append(min(allx) - margin)
    bbox.append(max(allx) + margin)
    bbox.append(min(ally) - margin)
    bbox.append(max(ally) + margin)
    return(bbox)
