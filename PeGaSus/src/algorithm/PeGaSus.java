package algorithm;


import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

public class PeGaSus{
    public long[] edges;
    public Int2IntOpenHashMap ID2IDX;
    public int numSubnodes;
    public int numSubedges;
    public int numSupernodes, numSuperedges;
    public Int2IntOpenHashMap subIDX2ID;
    public IntArrayList snList;
    public int[] rep;
    public double[] superRE;

    public IntArrayList[] insideSupernode;

    public double[] superPiW;
    public double[] superSquareW;
    public Int2DoubleOpenHashMap[] superEdgeCntW;

    public IntArrayList testNodeSetID;
    public IntArrayList testNodeSet;
    public IntArrayList targetNodeSetID;
    public IntArrayList targetNodeSet;


    public double targetRatio;

    public int divideMaxStep = 10;
    public int divThreshold = 500;
    public int T = 20;
    public int[] subHash;
    public int[] modelCost;
    public ObjectArrayList<IntArrayList> candidateDisjointGroups;

    int isolatedSupernode=-1;

    long[][] dup = new long[divThreshold][divThreshold];
    long cntFlag;

    PriorityQueue<Double> stayValue = new PriorityQueue<>(Collections.reverseOrder());

    public double alpha;
    public double[] alphaList;

    double log2(double a){
        return Math.log10(a)/Math.log10(2);
    }

    double edgeCost;
    public double getSuperedgeCost(int a){
        return 2*log2(numSupernodes-a);
    }


    public PeGaSus(IntArrayList _targetNodeSet, double _targetRatio, double _alpha, IntArrayList _testNodeSet){
        targetNodeSetID = _targetNodeSet.clone();
        if(_testNodeSet != null){
            testNodeSetID = _testNodeSet.clone();
        }
        else{
            testNodeSetID = _targetNodeSet.clone();
        }
        targetRatio = _targetRatio;
        cntFlag = 0;

        alpha =_alpha;
    }

    public double getInitialSize(){
        return numSubedges * 2 * log2(numSubnodes);
    }

    public double getSummarySize(){
        return numSuperedges * 2 * log2(numSupernodes) + numSubnodes * log2(numSupernodes);
    }

    public void addNode(int id){
        int idx = ID2IDX.getOrDefault(id, -1);
        if(idx < 0) {
            ID2IDX.put(id, numSubnodes);
            numSubnodes++;
        }
    }

    public boolean checkFormat(String number) {
        boolean isInt = false;
        try {
            Integer.parseInt(number);
            isInt = true;
        } catch (NumberFormatException e) {
        }
        return isInt;
    }


    public Int2IntOpenHashMap[] initializeSummary(String dataPath){
        String line;
        String[] parts;
        int cnt  = 0;
        ID2IDX = new Int2IntOpenHashMap();
        try {
            int srcID, dstID;
            BufferedReader br = new BufferedReader(new FileReader(dataPath));
            while ((line = br.readLine()) != null) {
                parts = line.split("\\s");
                if ((parts.length >= 2) && checkFormat(parts[0]) && checkFormat(parts[1])) {
                    try {
                        srcID = Integer.parseInt(parts[0]);
                        addNode(srcID);
                        dstID = Integer.parseInt(parts[1]);
                        addNode(dstID);
                        cnt++;
                    } catch (NumberFormatException e) {
                        System.out.println(e);
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            System.err.println(e);
        }

        targetNodeSet = new IntArrayList();
        testNodeSet = new IntArrayList();
        for(int t : targetNodeSetID) targetNodeSet.add(ID2IDX.get(t));
        for(int t : testNodeSetID) testNodeSet.add(ID2IDX.get(t));

        edges = new long[cnt];

        // Subgraph
        Int2IntOpenHashMap[] subAdj = new Int2IntOpenHashMap[numSubnodes];
        for (int i = 0; i < numSubnodes; i++) {
            subAdj[i] = new Int2IntOpenHashMap();
        }

        try {
            int srcIDX, dstIDX;
            BufferedReader br = new BufferedReader(new FileReader(dataPath));
            while ((line = br.readLine()) != null) {
                parts = line.split("\\s");
                if ((parts.length >= 2) && checkFormat(parts[0]) && checkFormat(parts[1])) {
                    try {
                        srcIDX = ID2IDX.get(Integer.parseInt(parts[0]));
                        dstIDX = ID2IDX.get(Integer.parseInt(parts[1]));
                        if (srcIDX == dstIDX) continue;
                        if (!subAdj[srcIDX].containsKey(dstIDX)) {
                            edges[numSubedges++] = (((long)srcIDX) << 32) + dstIDX;
                            subAdj[srcIDX].put(dstIDX, 1); subAdj[dstIDX].put(srcIDX, 1);
                        }
                    } catch (Exception e) {
                        System.out.println(e);
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            System.err.println(e);
        }

        System.out.println("|V|:\t" + numSubnodes);
        System.out.println("|E|:\t" + numSubedges);

        int id, idx;
        subIDX2ID = new Int2IntOpenHashMap();
        for(Int2IntOpenHashMap.Entry X : ID2IDX.int2IntEntrySet()){
            id = X.getIntKey();
            idx = X.getIntValue();
            subIDX2ID.put(idx, id);
        }

        alphaList = new double[1000];
        for(int i= 0;i<1000;i++) alphaList[i] = Math.pow(alpha, i);


        edgeCost = 2*log2(numSubnodes);

        snList = new IntArrayList(numSubnodes);
        superRE = new double[numSubnodes];
        modelCost = new int[numSubnodes];
        superPiW = new double[numSubnodes];
        superEdgeCntW = new Int2DoubleOpenHashMap[numSubnodes];
        insideSupernode = new IntArrayList[numSubnodes];
        superSquareW = new double[numSubnodes];
        rep = new int[numSubnodes];


        numSupernodes = numSubnodes;
        numSuperedges = numSubedges;

        for (int i = 0; i < numSubnodes; i++) {
            superRE[i] = 0;
            modelCost[i] = subAdj[i].keySet().size();
            rep[i] = i;
            snList.add(i);
            superEdgeCntW[i] = new Int2DoubleOpenHashMap();
            insideSupernode[i] = new IntArrayList(new int[]{i});
        }

        Int2IntOpenHashMap distance = new Int2IntOpenHashMap();
        Int2IntOpenHashMap hopHashMap = getHop(subAdj, targetNodeSet);

        for(int j = 0; j<numSubnodes; j++){
            distance.put(j, hopHashMap.get(j));
        }

        double sumPi = 0;
        double sumSqrPi = 0;
        for(int i =0;i<numSubnodes;i++) {
            int dis = distance.get(i);
            sumPi += 1./alphaList[dis];
            superPiW[i] = 1./alphaList[dis];
            superSquareW[i] = 1./alphaList[dis*2];
            sumSqrPi+=1./alphaList[dis*2];
        }
        double z = ((numSubnodes)/(sumPi*sumPi-sumSqrPi));
        z *= (numSubnodes-1);

        for(int i = 0;i<numSubnodes;i++){
            superPiW[i] *=  Math.sqrt(z);
            superSquareW[i] *=  z;
        }

        double eW;
        for(int i =0;i<numSubnodes;i++) {
            for(int nbd : subAdj[i].keySet()){
                eW = superPiW[i]*superPiW[nbd];
                superEdgeCntW[i].put(nbd, eW);
            }
        }


        return subAdj;
    }

    public Int2IntOpenHashMap getHop(Int2IntOpenHashMap[] subAdj, IntArrayList startIdxList){
        int maxHop=0;
        Int2IntOpenHashMap hop = new Int2IntOpenHashMap(0);
        for(int j = 0; j<numSubnodes; j++) hop.put(j, Integer.MAX_VALUE);

        Queue<Integer> q = new LinkedList<>();

        for(int s : startIdxList) hop.put(s, 0);

        q.addAll(startIdxList);

        while(!q.isEmpty()){
            int idx = q.poll();

            for(int nbd : subAdj[idx].keySet()){
                if(hop.get(nbd)==Integer.MAX_VALUE){
                    q.add(nbd);
                    hop.put(nbd,hop.get(idx)+1);
                    if(hop.get(idx)+1 > maxHop) maxHop=hop.get(idx)+1;
                }
            }
        }
        for(Int2IntOpenHashMap.Entry x : hop.int2IntEntrySet()){
            if(x.getIntValue()==Integer.MAX_VALUE){
                hop.put(x.getIntKey(), 2*maxHop);
            }
        }
        return hop;
    }
    public double getPiW(int A, int B){
        if(A==B) {
            double piW = superPiW[A]*superPiW[A]-superSquareW[A];
            piW = (piW>=0 ? piW: 0);
            return piW/2;
        }
        else return superPiW[A]*superPiW[B];
    }
    public double getSuperedgeCntW(int A, int B){
        double edgeCntW = superEdgeCntW[A].getOrDefault(B, 0);
        return Math.abs(edgeCntW);
    }

    public double getRE(int A, int B){
        double piW = getPiW(A,B);
        double edgeCntW = getSuperedgeCntW(A,B);

        return (superEdgeCntW[A].getOrDefault(B,0) > 0 ? edgeCost*(piW - edgeCntW) : edgeCost*edgeCntW);
    }

    public IntOpenHashSet getMergeCandidate(int A, int B){
        IntOpenHashSet candidate = new IntOpenHashSet();
        candidate.addAll(superEdgeCntW[A].keySet());
        candidate.addAll(superEdgeCntW[B].keySet());

        if(candidate.contains(B)) {
            candidate.add(A);
            candidate.remove(B);
        }

        return candidate;
    }

    public double getMergeCost(int A, int B){
        IntOpenHashSet candidate = getMergeCandidate(A,B);
        double edgeCntW, piW, denseCost, sparseCost;
        double mergeCost=0;

        if(candidate.contains(A)){
            edgeCntW = getSuperedgeCntW(A,A) + getSuperedgeCntW(A,B) + getSuperedgeCntW(B,B);
            piW = Math.pow(superPiW[A]+superPiW[B],2) - (superSquareW[A]+superSquareW[B]);
            piW = (piW>=0 ? piW/2 : 0);

            denseCost = edgeCost*(piW-edgeCntW) + getSuperedgeCost(1);


            sparseCost = edgeCost*edgeCntW;

            mergeCost += Math.min(denseCost, sparseCost);
        }

        for(int can : candidate){
            if(can == A) continue;
            edgeCntW = getSuperedgeCntW(A,can) + getSuperedgeCntW(B, can);
            piW = getPiW(A,can) + getPiW(B, can);
            denseCost = edgeCost*(piW-edgeCntW) + getSuperedgeCost(1);
            sparseCost = edgeCost*edgeCntW;

            mergeCost += Math.min(denseCost, sparseCost);
        }

        return mergeCost;
    }

    public long getBestPair(int[] sunList, double threshold,int length){
        long uniqValue = cntFlag;
        int sz = length;
        if(sz < 2) return -1;

        Random rand = new Random();
        int randomN;
        long SRCDST;
        double bestValue;
        long maxReductionPair;

        bestValue = -Double.MAX_VALUE;
        maxReductionPair=-1;

        int A,B;
        int idxA, idxB;
        for(int i =0;i<sz;i++){
            randomN = rand.nextInt(sz*(sz-1));
            A = randomN / (sz-1); B = randomN % (sz-1);
            if(A<=B) B++;
            if(dup[A][B] == uniqValue) continue;
            dup[A][B] = uniqValue;

            idxA = sunList[A]; idxB = sunList[B];
            SRCDST = ((long)A << 32) + B;

            double beforeRE = superRE[idxA]+superRE[idxB]-getRE(idxA,idxB);
            beforeRE += (modelCost[idxA]+modelCost[idxB]-(superEdgeCntW[idxA].getOrDefault(idxB,0)>0 ? 1: 0))*getSuperedgeCost(0);

            double afterRE = getMergeCost(idxA, idxB);

            double transError = 1 -(afterRE/beforeRE);

            if(bestValue < transError) {
                maxReductionPair = SRCDST;
                bestValue =  transError;
            }

        }

        if(bestValue >= threshold){
            return maxReductionPair;
        }
        else{
            stayValue.add(bestValue);
            return -1;
        }
    }

    public int[][] signature(int k){
        final int maxCore = Runtime.getRuntime().availableProcessors();
        ExecutorService executorService = Executors.newFixedThreadPool(maxCore);
        List<Future<int[]>> resultList = new ArrayList<>();
        int[][] ans = new int[numSupernodes][k+1];


        for (int i=0; i<k; i++)
        {
            MinH calculate = new MinH(numSupernodes, numSubnodes, snList, insideSupernode, edges);
            Future<int[]> result = executorService.submit(calculate);
            resultList.add(result);
        }


        executorService.shutdown();
        try {
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
        }

        int col = 0;
        for(Future<int[]> future : resultList)
        {
            try
            {
                int j = 0;
                for(int i:future.get()){
                    ans[j++][col] = i;
                }
            }
            catch (InterruptedException | ExecutionException e)
            {
                e.printStackTrace();
            }
            col++;
        }

        for(int i=0;i<numSupernodes;i++) {
            ans[i][col] = snList.getInt(i);
        }

        int [][] sorted = Arrays.stream(ans).sorted(Arrays::compare).toArray(int[][]::new);

        return sorted;
    }

    // Divide signature matrix into supernode groups (spec. 500)
    public ArrayList<int[]> divSupernode(int[][] array, int k){
        ArrayList<int[]> tmp;


        int[] supernode = new int[numSupernodes];


        for(int i=0;i<numSupernodes;i++){
            supernode[i] = array[i][0];
        }

        tmp = listN(supernode, 0);

        return recursiveDiv(array, tmp, 0, k);
    }

    public int[] column(int[][] array, int col, int first, int last){
        int[] ans = new int[last-first+1];
        int index = 0;
        for(int i = first; i<=last; i++){
            ans[index++] = array[i][col];
        }
        return ans;
    }

    public ArrayList<int[]> recursiveDiv(int[][] array, ArrayList<int[]> tmp, int col, int k){
        ArrayList<int[]> ans = new ArrayList<>();

        for(int[] x:tmp){
            int first = x[0];
            int last = x[1];
            int length = last - first + 1;
            if((length>divThreshold)&&(col<k)){
                col++;
                int [] tempp = column(array, col, first, last);

                ArrayList<int[]> tempo = listN(tempp, first);

                ans.addAll(recursiveDiv(array, tempo ,col, k));
                col--;
            }

            // Groups are not smaller than divThreshold
            if((length>divThreshold)&&(col==k)){
                // Divide each group by using first & Last
                int quot = length / divThreshold;
                int tmpLast = first + divThreshold -1;
                quot--;
                while(tmpLast<last){
                    ans.add(column(array, k, first, tmpLast));
                    first += divThreshold;
                    if(quot!=0){
                        tmpLast += divThreshold;
                        quot--;
                        continue;
                    }
                    tmpLast = last;
                }
                // Left over
                ans.add(column(array, k, first, tmpLast));
            }
            if(length<=divThreshold) {
                ans.add(column(array, k, first, last));
            }
        }
        return ans;
    }

    public ArrayList<int[]> listN(int[] array, int first){
        ArrayList<int[]> ans = new ArrayList<>();
        int tmp = array[0];
        int last = first-1;

        for(int i=0;i<array.length;i++){
            if(tmp==array[i]) last++;
            else{
                int[] model = new int[2];
                model[0] = first;
                model[1] = last;
                ans.add(model);
                tmp=array[i];
                first = last++ +1;
            }
        }
        int[] model = new int[2];
        model[0] = first;
        model[1] = last;
        ans.add(model);
        return ans;
    }

    public void merge(ArrayList<int[]> candidateSNGroup, double threshold) {
        for (int[] SL : candidateSNGroup) {
            int szSupernodeSet = SL.length;

            if(szSupernodeSet <2) continue;

            int numSkips = 0;
            while(numSkips < Math.max(log2(szSupernodeSet), 1)){
                cntFlag += 1;

                long bestPair = getBestPair(SL, threshold, szSupernodeSet);
                if (bestPair == -1) {
                    numSkips += 1;
                }
                else{
                    int A = (int) (bestPair >> 32); int B = (int) (bestPair & 0x7FFFFFFFL);
                    int idxA = SL[A]; int idxB = SL[B];

                    mergeInner(idxA, idxB);
                    szSupernodeSet -= 1;
                    SL[B] = SL[szSupernodeSet];
                    numSkips = 0;
                    double SR = getSummarySize()/getInitialSize();


                    if (targetRatio > (getSummarySize() / getInitialSize())) {
                        break;
                    }
                }
                if (targetRatio > (getSummarySize() / getInitialSize())) {
                    break;
                }
            }
        }
    }

    public void summarize(){
        subHash = new int[numSubnodes];
        candidateDisjointGroups = new ObjectArrayList<>(0);

        double threshold;
        double it = 1 ;

        int thresholdIdx;
        double thresholdRatio =  0.1;
        threshold = 0.5;
        while(true){
//            System.out.println("Iter\t"+it+"\tCompression Ratio\t"+(getSummarySize()/getInitialSize()));
            int[][] sorted = signature(divideMaxStep);
            ArrayList<int[]> tmp = divSupernode(sorted, divideMaxStep);
            merge(tmp, threshold);

            if (targetRatio > (getSummarySize() / getInitialSize())) {
                stayValue = null;
                break;
            }

            if(stayValue.size()>1){
                thresholdIdx = (int)(stayValue.size() * thresholdRatio);
                for(int i = 0 ; i<thresholdIdx;i++) stayValue.remove();
                threshold = stayValue.remove();
            }
            if(threshold<0) threshold= 0;

            stayValue.clear();

            if(it == T) {
                stayValue = null;
                break;
            }
            it++;
        }

        for(int n: snList){
            superEdgeCntW[n].int2DoubleEntrySet().removeIf(e -> (e.getDoubleValue() < 0));
        }

        int i = 0;
        while(i < snList.size()){
            int sn = snList.getInt(i);
            i++;
            if(superEdgeCntW[sn].keySet().size() == 0){
                if(isolatedSupernode == -1){
                    modelCost[sn] = 0;
                    isolatedSupernode =  sn;
                }
                else{
                    insideSupernode[isolatedSupernode].addAll(insideSupernode[sn]);
                    superEdgeCntW[sn]=null;
                    insideSupernode[sn]=null;
                    modelCost[sn] = -Integer.MAX_VALUE;

                    numSupernodes--;
                    snList.set(rep[sn], snList.getInt(numSupernodes));
                    rep[snList.getInt(rep[sn])] = rep[sn];
                    rep[sn] = -1;
                    snList.popInt();
                    i--;
                }
            }
        }

        if (targetRatio > (getSummarySize() / getInitialSize())) {
            return;
        }

        if (targetRatio < (getSummarySize() / getInitialSize())) {
            System.out.println("TopNDrop");
            topNDrop();
        }


    }

    public double getDropError(int A, int B){
        double piW = getPiW(A,B);
        double edgeCntW = getSuperedgeCntW(A,B);

        return 2*edgeCntW - piW;
    }

    public double topNDrop(){
        double evalConsumeTime = 0;
        int SRC, DST;

        int[][] sedges = new int[numSuperedges][2];
        double[] cost = new double[numSuperedges];
        int[] idxs = new int[numSuperedges];
        int[] medians = new int[numSuperedges];

        int col = 0;

        for(int n: snList) {
            for (Int2DoubleOpenHashMap.Entry e : Int2DoubleMaps.fastIterable(superEdgeCntW[n])) {
                if (e.getIntKey() >= n) {
                    cost[col] = getDropError(e.getIntKey(), n);
                    sedges[col][0] = n;
                    sedges[col][1] = e.getIntKey();
                    idxs[col] = col;
                    col++;
                }
            }
        }
        double summarySize = getSummarySize();
        int ss = 0, ee = numSuperedges-1;

        while(summarySize > targetRatio*getInitialSize()) {
            int edgeLft;
            edgeLft = numSuperedges - (int)((targetRatio*getInitialSize() - numSubnodes * log2(numSupernodes)) / (2 * log2(numSupernodes)));

            int _tmp;
            for (int i = ss; i <= ee; i++) {
                medians[i - ss] = idxs[i];
            }
            int chosen = -1;
            int mediansLeft = ee - ss + 1;
            while (mediansLeft > 1) {
                for (int j = 0; j < mediansLeft; j += 5) {
                    int mlen = (j + 5 < mediansLeft) ? 5 : (mediansLeft - j);
                    if (mlen == 5) {
                        if (cost[medians[j]] > cost[medians[j + 1]]) {
                            _tmp = medians[j];
                            medians[j] = medians[j + 1];
                            medians[j + 1] = _tmp;
                        }
                        if (cost[medians[j + 2]] > cost[medians[j + 3]]) {
                            _tmp = medians[j + 2];
                            medians[j + 2] = medians[j + 3];
                            medians[j + 3] = _tmp;
                        }
                        if (cost[medians[j + 1]] < cost[medians[j + 3]]) {
                            medians[j + 3] = medians[j + 4];
                            if (cost[medians[j + 2]] > cost[medians[j + 3]]) {
                                _tmp = medians[j + 2];
                                medians[j + 2] = medians[j + 3];
                                medians[j + 3] = _tmp;
                            }
                        } else {
                            medians[j + 1] = medians[j + 4];
                            if (cost[medians[j]] > cost[medians[j + 1]]) {
                                _tmp = medians[j];
                                medians[j] = medians[j + 1];
                                medians[j + 1] = _tmp;
                            }
                        }
                        if (cost[medians[j + 1]] < cost[medians[j + 3]]) {
                            chosen = (cost[medians[j + 1]] > cost[medians[j + 2]]) ? medians[j + 1] : medians[j + 2];
                        } else {
                            chosen = (cost[medians[j]] > cost[medians[j + 3]]) ? medians[j] : medians[j + 3];
                        }
                    } else if (mlen == 4) {
                        if (cost[medians[j]] > cost[medians[j + 1]]) {
                            _tmp = medians[j];
                            medians[j] = medians[j + 1];
                            medians[j + 1] = _tmp;
                        }
                        if (cost[medians[j + 2]] > cost[medians[j + 3]]) {
                            _tmp = medians[j + 2];
                            medians[j + 2] = medians[j + 3];
                            medians[j + 3] = _tmp;
                        }
                        if (cost[medians[j]] < cost[medians[j + 2]]) {
                            chosen = (cost[medians[j + 1]] < cost[medians[j + 2]]) ? medians[j + 1] : medians[j + 2];
                        } else {
                            chosen = (cost[medians[j]] > cost[medians[j + 3]]) ? medians[j] : medians[j + 3];
                        }
                    } else if (mlen == 3) {
                        if (cost[medians[j]] > cost[medians[j + 1]]) {
                            _tmp = medians[j];
                            medians[j] = medians[j + 1];
                            medians[j + 1] = _tmp;
                        }
                        if (cost[medians[j]] < cost[medians[j + 2]]) {
                            chosen = (cost[medians[j + 1]] < cost[medians[j + 2]]) ? medians[j + 1] : medians[j + 2];
                        } else chosen = medians[j];
                    } else if (mlen == 2) {
                        chosen = (cost[medians[j]] < cost[medians[j + 1]]) ? medians[j] : medians[j + 1];
                    } else {
                        chosen = medians[j];
                    }
                    medians[j / 5] = chosen;
                }
                mediansLeft = (mediansLeft + 4) / 5;

            }
            chosen = medians[0];
            int _ss = ss, _ee = ee, _mid = ee;
            while (_ss <= _mid) {
                if (cost[idxs[_ss]] < cost[chosen]) {
                    _ss++;
                } else if (cost[idxs[_ss]] > cost[chosen]) {
                    _tmp = idxs[_ss];
                    idxs[_ss] = idxs[_mid];
                    idxs[_mid] = idxs[_ee];
                    idxs[_ee] = _tmp;
                    _ee--;
                    _mid--;
                } else {
                    _tmp = idxs[_ss];
                    idxs[_ss] = idxs[_mid];
                    idxs[_mid] = _tmp;
                    _mid--;
                }
            }
            if (edgeLft <= (_ss - ss)) {
                ee = _ss - 1;
            } else if (edgeLft <= (_ee - ss + 1)) {
                for (int i = ss; i < ss + edgeLft; i++) {
                    SRC = sedges[idxs[i]][0];
                    DST = sedges[idxs[i]][1];

                    numSuperedges--;
                    modelCost[SRC]--;
                    if(SRC!=DST) modelCost[DST]--;
                    superEdgeCntW[SRC].remove(DST);
                    superEdgeCntW[DST].remove(SRC);

                    if(modelCost[SRC] == 0){
                        if(isolatedSupernode==-1){
                            isolatedSupernode = SRC;
                        }
                        else{
                            insideSupernode[isolatedSupernode].addAll(insideSupernode[SRC]);
                            superEdgeCntW[SRC]=null;
                            insideSupernode[SRC]=null;
                            modelCost[SRC] = -Integer.MAX_VALUE;

                            numSupernodes--;
                            snList.set(rep[SRC], snList.getInt(numSupernodes));
                            rep[snList.getInt(rep[SRC])] = rep[SRC];
                            rep[SRC] = -1;
                            snList.popInt();
                        }
                    }
                    if(SRC != DST && modelCost[DST] == 0){
                        if(isolatedSupernode==-1){
                            isolatedSupernode = DST;
                        }
                        else{
                            insideSupernode[isolatedSupernode].addAll(insideSupernode[DST]);
                            superEdgeCntW[DST]=null;
                            insideSupernode[DST]=null;
                            modelCost[DST] = -Integer.MAX_VALUE;

                            numSupernodes--;
                            snList.set(rep[DST], snList.getInt(numSupernodes));
                            rep[snList.getInt(rep[DST])] = rep[DST];
                            rep[DST] = -1;
                            snList.popInt();
                        }
                    }

                    if (targetRatio > (getSummarySize() / getInitialSize())) {
                        break;
                    }
                }
                break;
            } else {
                for (int i = ss; i <= _ee; i++) {
                    SRC = sedges[idxs[i]][0];
                    DST = sedges[idxs[i]][1];

                    numSuperedges--;
                    modelCost[SRC]--;
                    if (SRC != DST) modelCost[DST]--;
                    superEdgeCntW[SRC].remove(DST);
                    superEdgeCntW[DST].remove(SRC);

                    if (modelCost[SRC] == 0) {
                        if (isolatedSupernode == -1) {
                            isolatedSupernode = SRC;
                        } else {
                            insideSupernode[isolatedSupernode].addAll(insideSupernode[SRC]);
                            superEdgeCntW[SRC] = null;
                            insideSupernode[SRC] = null;
                            modelCost[SRC] = -Integer.MAX_VALUE;

                            numSupernodes--;
                            snList.set(rep[SRC], snList.getInt(numSupernodes));
                            rep[snList.getInt(rep[SRC])] = rep[SRC];
                            rep[SRC] = -1;
                            snList.popInt();
                        }
                    }
                    if (SRC != DST && modelCost[DST] == 0) {
                        if (isolatedSupernode == -1) {
                            isolatedSupernode = DST;
                        } else {
                            insideSupernode[isolatedSupernode].addAll(insideSupernode[DST]);
                            superEdgeCntW[DST] = null;
                            insideSupernode[DST] = null;
                            modelCost[DST] = -Integer.MAX_VALUE;

                            numSupernodes--;
                            snList.set(rep[DST], snList.getInt(numSupernodes));
                            rep[snList.getInt(rep[DST])] = rep[DST];
                            rep[DST] = -1;
                            snList.popInt();
                        }
                    }

                    if (targetRatio > (getSummarySize() / getInitialSize())) {
                        break;
                    }
                }
                if (targetRatio > (getSummarySize() / getInitialSize())) {
                    break;
                }
                ss = _ee + 1;
            }
        }
        return evalConsumeTime;
    }

    public void editSuperedge(int A, int B, boolean add, double edgeCntW){
        if(add){
            if(superEdgeCntW[A].getOrDefault(B,0) <= 0) {
                numSuperedges++;
                modelCost[A]++;
                if(A!=B) modelCost[B]++;

            }

            superEdgeCntW[A].put(B, edgeCntW);
            superEdgeCntW[B].put(A, edgeCntW);
        }
        else{
            if(superEdgeCntW[A].getOrDefault(B,0) > 0) {
                numSuperedges--;
                modelCost[A]--;
                if(A!=B) modelCost[B]--;
            }

            superEdgeCntW[A].put(B, -edgeCntW);
            superEdgeCntW[B].put(A, -edgeCntW);
        }
    }


    public void mergeInner(int A, int B){
        IntOpenHashSet candidate = getMergeCandidate(A,B);

        double piW, edgeCntW, denseCost, sparseCost, reAPS = 0;
        boolean dense;

        if(candidate.contains(A)){
            edgeCntW = getSuperedgeCntW(A,A) + getSuperedgeCntW(A,B) + getSuperedgeCntW(B,B);
            piW = Math.pow(superPiW[A]+superPiW[B],2) - (superSquareW[A]+superSquareW[B]);
            piW = (piW>=0 ? piW/2 : 0);
            denseCost = edgeCost*(piW-edgeCntW);
            sparseCost = edgeCost*edgeCntW;


            dense = (denseCost+ getSuperedgeCost(1)) < sparseCost;

            if(dense){
                reAPS += denseCost;
            }
            else{
                reAPS += sparseCost;
            }

            editSuperedge(A,A, dense, edgeCntW);

            if(superEdgeCntW[A].get(B) > 0) {
                numSuperedges--;
                modelCost[A]--;
                if(A!=B) modelCost[B]--;
            }
            superEdgeCntW[A].remove(B);
        }

        for(int can : candidate){
            if(can == A) continue;
            superRE[can] -= getRE(A,can);
            superRE[can] -= getRE(B,can);
            piW = getPiW(A,can) + getPiW(B, can);
            edgeCntW = getSuperedgeCntW(A, can) + getSuperedgeCntW(B, can);

            denseCost = edgeCost*(piW-edgeCntW);
            sparseCost = edgeCost*edgeCntW;


            dense = (denseCost+ getSuperedgeCost(1)) < sparseCost;

            if(dense){
                reAPS += denseCost;
                superRE[can] += denseCost;
            }
            else{
                reAPS += sparseCost;
                superRE[can] += sparseCost;
            }

            editSuperedge(A,can, dense, edgeCntW);

            if(superEdgeCntW[B].get(can) > 0) {
                numSuperedges--;
                modelCost[B]--;
                if(B!=can) modelCost[can]--;
            }
            superEdgeCntW[can].remove(B);
        }

        if(superEdgeCntW[B].get(B) > 0) {
            numSuperedges--;
            modelCost[B]--;
        }

        insideSupernode[A].addAll(insideSupernode[B]);
        superPiW[A] += superPiW[B];
        superSquareW[A] += superSquareW[B];

        superEdgeCntW[B] = null;

        insideSupernode[B]=null;

        superPiW[B] = 0;
        superSquareW[B]=0;

        numSupernodes--;
        snList.set(rep[B], snList.getInt(numSupernodes));
        rep[snList.getInt(rep[B])] = rep[B];
        rep[B] = -1;
        snList.popInt();

        superRE[A] = reAPS; superRE[B] = -Double.MAX_VALUE;
        modelCost[B] = -Integer.MAX_VALUE;
    }

}
