import algorithm.PeGaSus;
import it.unimi.dsi.fastutil.ints.*;

import java.io.*;
import java.util.LinkedList;
import java.util.Queue;
import java.util.stream.IntStream;

public class RUN {
    public static void main(String[] args) {
        final String dataPath = args[0];
        System.out.println("dataPath: " + dataPath);

        double targetRatio = Double.parseDouble(args[1]);
        System.out.println("Target Ratio: " + targetRatio);

        String targetNodeSetFileName = args[2];
        System.out.println("Target Node Set File Name: " + targetNodeSetFileName);

        boolean isCheckPersonalizedErr = Boolean.parseBoolean(args[3]);
        System.out.println("is Checking Personalized Error: " + isCheckPersonalizedErr);

        boolean isSaveSummaryGraph = Boolean.parseBoolean(args[4]);
        System.out.println("is Saving Summary Graph: " + isSaveSummaryGraph);

        double alpha = Double.parseDouble(args[5]);
        System.out.println("Alpha: " + alpha);

        double start = System.currentTimeMillis();
        IntArrayList targetNodeSet = readSeedNodeSet(targetNodeSetFileName);
        if(targetNodeSet.size() == 0) {
            System.out.println("Check Target Node Set File");
            return;
        }

        // If you want to compute a personalized error for some node, you replace null to IntArrayList containing the nodes.
        PeGaSus m = new PeGaSus(targetNodeSet, targetRatio, alpha, null);

        Int2IntOpenHashMap[] subAdj = null;
        if(isCheckPersonalizedErr) {
            subAdj = m.initializeSummary(dataPath);
        }
        else{
            m.initializeSummary(dataPath);
        }
        m.summarize();
        System.out.println("Compression Ratio\t"+targetRatio+"\tExecution time\t " + (System.currentTimeMillis() - start) / 1000.0 + "s");
        if(isCheckPersonalizedErr) System.out.println("Personalized Error\t"+streamPersonalizedErrorforSeedNodeSet(m, subAdj, m.targetNodeSet));
        if(isSaveSummaryGraph) saveSummaryGraph(m);
    }

    static public void saveSummaryGraph(PeGaSus sumG) {
        String outputPath = "./output/";
        File f = new File(outputPath);

        if (!f.exists()) {
            try{
                f.mkdirs();
            }
            catch(Exception e){
                e.getStackTrace();
            }
        }
        f = new File(outputPath + File.separator+"summary_graph.txt");
        try {
            FileWriter fw = new FileWriter(f);
            fw.write("<Subnode of each supernode>");
            fw.write(System.getProperty( "line.separator" ));
            for (int sup_v : sumG.snList) {
                fw.write(String.format("%d", sup_v));
                for (int subIdx : sumG.insideSupernode[sup_v]) {
                    fw.write("\t"+String.format("%d", sumG.subIDX2ID.get(subIdx)));
                }
                fw.write(System.getProperty( "line.separator" ));
            }
            fw.write("<Superedge info>");
            fw.write(System.getProperty( "line.separator" ));
            for(int sup_v: sumG.snList){
                for(Int2DoubleOpenHashMap.Entry x : sumG.superEdgeCntW[sup_v].int2DoubleEntrySet()){
                    int NBD = x.getIntKey();
                    double edgeCntW = x.getDoubleValue();
                    if(edgeCntW>0){
                        if(sup_v>=NBD){
                            fw.write(String.format("%d", sup_v) +"\t"+ String.format("%d", NBD));
                            fw.write(System.getProperty( "line.separator" ));
                        }
                    }
                }
            }
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static public double streamPersonalizedErrorforSeedNodeSet(PeGaSus sumG, Int2IntOpenHashMap[] subAdj, IntArrayList testNodeSet){
        double PE=0;
        IntStream stream;

        stream = IntStream.range(0, testNodeSet.size());
        // If you want to compute a personalized loss fast, use .parallel()
//        PE += stream.parallel().mapToDouble(idx -> _personalizedErrorforSeedNodeSet(subAdj, idx)).sum();
        PE += stream.mapToDouble(idx -> _personalizedErrorforSeedNodeSet(sumG, subAdj, testNodeSet.getInt(idx))).sum();
        PE /= testNodeSet.size();

        return PE;
    }

    public static Int2IntOpenHashMap getHop(Int2IntOpenHashMap[] subAdj, int startIdx){
        int maxHop =0;
        Int2IntOpenHashMap hop = new Int2IntOpenHashMap(0);
        for(int j = 0; j<subAdj.length; j++) hop.put(j, Integer.MAX_VALUE);

        Queue<Integer> q = new LinkedList<>();
        hop.put(startIdx, 0);

        q.add(startIdx);

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

    public static double _personalizedErrorforSeedNodeSet(PeGaSus sumG, Int2IntOpenHashMap[] subAdj, int u){
        double error;
        double semiErr = 0;

        Int2IntOpenHashMap snsHop;

        int numSubnodes = subAdj.length;

        double maxError;
        Int2IntOpenHashMap hopMap = getHop(subAdj, u);
        double piA = 0; double piB=0;
        for(int idx=0;idx<numSubnodes;idx++) {
            int h = hopMap.get(idx);
            piA += 1./sumG.alphaList[h];
            piB += 1./sumG.alphaList[2*h];
        }
        maxError = (piA*piA - piB)/2;

        int temp;
        int[] inv = new int[numSubnodes];
        for(int v: sumG.snList){
            for(int _u: sumG.insideSupernode[v]){
                inv[_u] = v;
            }
        }

        snsHop = getHop(subAdj, u);
        error = 0.;

        for(long e : sumG.edges){
            int src = (int)(e >> 32);
            int dst = (int)(e & 0x7FFFFFFFL);

            int SRC = inv[src];
            int DST = inv[dst];

            if(sumG.superEdgeCntW[SRC].getOrDefault(DST,0)<=0) error += ((1/(sumG.alphaList[snsHop.get(src)+snsHop.get(dst)])));
        }
        semiErr = 0;
        for(int A: sumG.snList){
            for(Int2DoubleOpenHashMap.Entry X:sumG.superEdgeCntW[A].int2DoubleEntrySet()){
                int B = X.getIntKey();
                if (A>B) continue;
                if(X.getDoubleValue()>0){
                    for(int a : sumG.insideSupernode[A]){
                        temp = snsHop.get(a);
                        for(int b : sumG.insideSupernode[B]){
                            if(a!=b && subAdj[a].getOrDefault(b,0) == 0){
                                semiErr += ((1/(sumG.alphaList[temp+snsHop.get(b)])));
                            }
                        }
                    }
                    if(A==B) semiErr/=2;
                    error += semiErr;
                }
                semiErr=0;
            }
        }

        return error/maxError;
    }


    public static IntArrayList readSeedNodeSet(String dataPath){
        IntArrayList SNS;
        String line;

        SNS = new IntArrayList();

        try {
            BufferedReader br = new BufferedReader(new FileReader(dataPath));
            while ((line = br.readLine()) != null) {
                SNS.add(Integer.parseInt(line));
            }
            br.close();
        } catch (IOException e) {
            System.err.println(e);
        }



        return SNS;

    }


}
