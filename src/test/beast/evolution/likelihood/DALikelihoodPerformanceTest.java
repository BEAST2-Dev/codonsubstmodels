package test.beast.evolution.likelihood;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import beast.likelihood.DACodonTreeLikelihood;
import beast.tree.InternalNodeStates;
import beast.util.XMLParserException;
import codonmodels.CodonFrequencies;
import org.junit.Test;
import test.beast.evolution.CodonData;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Random;


/**
 * Performance tests between standard TreeLikelihood and
 * DACodonTreeLikelihood.
 *
 * @author Walter Xie
 */
public class DALikelihoodPerformanceTest {
    CodonAlignment codonAlignment;
    SiteModel siteModel;
    Tree tree;

    DecimalFormat df = new DecimalFormat("#.00");



    @Test
    public void benchmarkingCodon100(){
//        initF3X4(128, 100, true);
//        benchmarking(1);
        taxa2to128(100, 1);
    }

    /**
     * 10 iteration(s) :
     * test 2 taxa : DA likelihood 2.50 times faster
     * test 4 taxa : DA likelihood 4.43 times faster
     * test 8 taxa symmetric : DA likelihood 9.77 times faster
     * test 8 taxa asymmetric : DA likelihood 15.57 times faster
     * test 16 taxa : DA likelihood 19.92 times faster
     * test 32 taxa : DA likelihood 24.45 times faster
     */
    @Test
    public void benchmarkingCodon500(){
        taxa2to128(500, 1);
    }

    protected void taxa2to128(int nCodon, int iteration) {
        // speed symmetric < asymmetric
        String[] test = new String[]{"2 taxa", "4 taxa symmetric", "4 taxa asymmetric",
                "8 taxa symmetric", "8 taxa asymmetric", "16 taxa", "32 taxa", "64 taxa", "128 taxa"};
        boolean[] symmetric = new boolean[]{true, true, false, true, false, true, true, true, true};
        int[] nTaxa = new int[]{2, 4, 4, 8, 8, 16, 32, 64, 128};
        double[] faster = new double[test.length];
//        int ite = 10;

        for (int i = 0; i < test.length; i++) {
            initF3X4(nTaxa[i], nCodon, symmetric[i]);
            faster[i] = benchmarking(iteration);
        }

        System.out.println("\n=============== Summary " + nCodon + " codons ===============");
        System.out.println(iteration + " iteration(s) :\n");
        for (int i = 0; i < faster.length; i++) {
            System.out.println("test " + test[i] + " : DA likelihood " + df.format(faster[i]) + " times faster");
        }
    }

    // seq length = nCodon * 3
    private void initF3X4(int nTaxa, int nCodon, boolean symmetric) {
        codonAlignment = null;
        siteModel = null;
        tree = null;
        Alignment data = null;
        try {
            data = CodonData.getAlignment(nTaxa, nCodon);
        } catch (XMLParserException e) {
            e.printStackTrace();
        }

        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial", "verbose", false);

        String pi = "F3X4";
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", false);
        siteModel = CodonData.getSiteModel("0.3", "5", codonFreq, false);

        String newickTree;
        if (nTaxa == 2) { // NullPointerException in TreeLikelihood.setStates
            newickTree = "(t1:0.5,t2:0.1):0.0;";
        } else if (nTaxa == 4 && symmetric) {
            newickTree = "((t1:0.5,t2:0.1):0.2, (t3:0.5,t4:0.1):0.3);";
        } else if (nTaxa == 4) { // not symmetric
            newickTree = "(((t1:0.5,t2:0.1):0.2, t3:0.5):0.1,t4:0.1);";
        } else if (nTaxa == 8 && symmetric) {
            newickTree = "(((t1:0.5,t2:0.1):0.2, (t3:0.5,t4:0.1):0.3):0.1, ((t5:0.4,t6:0.2):0.3, (t7:0.3,t8:0.2):0.1):0.2);";
        } else if (nTaxa == 8) { // not symmetric
            newickTree = "(((((((t1:0.5,t2:0.1):0.2, t3:0.5):0.1,t4:0.1):0.2, t5:0.4):0.1,t6:0.2):0.1,t7:0.3):0.1,t8:0.2 );";
        } else if (nTaxa == 16) {
            newickTree = "((t6:0.8328606829,t5:0.07987240003):0.8731191335,(((((t10:0.4489309546,t2:0.5118521685):0.6122735683,t12:0.3223016059):0.7295662321,(t15:0.6618216026,(t4:0.6256095201,t1:0.7471578235):0.8436422099):0.3587985209):0.07577025215,(((t11:0.375966636,(t14:0.01772094774,t9:0.5878796382):0.5867906068):0.2359019916,t8:0.2657493246):0.675317778,t16:0.8670396307):0.3450425724):0.1524739896,(t3:0.4887036914,(t13:0.4686036934,t7:0.1729974856):0.1101265033):0.1990886456):0.8794986582);";
        } else if (nTaxa == 32) {
            newickTree = "(((t15:0.9542181483,t19:0.3497522059):0.3365416529,t2:0.5435624667):0.9214519702,((((((t21:0.7185004468,t20:0.5407348312):0.7264935684,(t17:0.5117207766,(t8:0.3790637162,(t30:0.7958646573,(t14:0.4777974519,t18:0.9716053747):0.7753737369):0.6706928122):0.6051157189):0.890223223):0.1895852804,(t22:0.4900219755,t13:0.05010107416):0.04109163722):0.06207286357,((t6:0.5607458302,t16:0.9570369916):0.4422796173,t31:0.8856992035):0.4335002408):0.7788049644,(((t7:0.7984650929,t4:0.9124820076):0.224632408,(t1:0.544289869,t32:0.702123753):0.5266149947):0.8062962941,((t24:0.9591187276,t27:0.7903348505):0.9343056851,t10:0.924266933):0.8483990482):0.8787279353):0.4969263494,(((t5:0.596928519,t23:0.6392580203):0.3214438602,t9:0.1940213139):0.2096933059,((t11:0.8579683306,(t28:0.9736292204,t26:0.06027402054):0.07778765121):0.1267342451,((t29:0.8047343672,t25:0.236599003):0.0191930905,(t12:0.9822205401,t3:0.1412437484):0.1372304766):0.6479954512):0.6275273843):0.5978206636):0.5197886776);";
        } else if (nTaxa == 64) {
            newickTree = "(((((((t48:0.1194739025,(t23:0.5121433961,t55:0.9955580584):0.8700787935):0.9662126682,(((t49:0.4517673494,((t56:0.451374223,t30:0.1933717055):0.9033906148,t33:0.02580949571):0.5926264862):0.5438068295,((t53:0.6748554085,(t11:0.05578524782,t9:0.06882344):0.6639610962):0.6540090321,(t20:0.3769058259,t24:0.4544169321):0.8009162853):0.5525573976):0.5693003454,(t28:0.4067790734,t27:0.7071143268):0.2404345684):0.7761483749):0.1718712663,(t41:0.7353895123,(t52:0.08006391395,t42:0.2578710266):0.2359069309):0.05422290228):0.1061627998,((t34:0.5635265361,t22:0.9561055887):0.09368016128,(t63:0.6213681796,(t14:0.02218521596,(t45:0.2115799016,(t10:0.8981468673,(t18:0.9982340881,t32:0.9741900661):0.459091333):0.09610739653):0.08724628296):0.371052):0.1051064976):0.5608417923):0.1209174024,(((((t36:0.5428988964,t44:0.335287529):0.9125117145,((t4:0.7913779092,t1:0.7296556616):0.6592618392,((t12:0.9368871376,t51:0.1481878199):0.05215751473,t21:0.6213180579):0.4304548244):0.5414799473):0.9699665923,((((((t29:0.7395233617,t25:0.9559190974):0.06876781909,t13:0.5378922906):0.8743523394,(t64:0.3304065655,(t5:0.9110768361,t40:0.1287771119):0.904186056):0.3743544784):0.3717798328,(((t60:0.1011154889,(t61:0.01538130431,(t50:0.1364783351,t62:0.05287580518):0.5683464571):0.9057637458):0.8347899029,(t15:0.8597799973,t3:0.3062982489):0.7550411909):0.5969822519,t35:0.164621867):0.1106226472):0.5239213563,t38:0.9545394897):0.01059990446,((t43:0.4596710789,t46:0.5547174369):0.9316371751,t6:0.6194477731):0.9156797142):0.6047292564):0.3184891245,((t31:0.1859473379,t7:0.3811759248):0.1085469136,(((t54:0.4991254942,t2:0.9472117536):0.793728224,t39:0.6302277825):0.1995151127,t19:0.8844033491):0.07947330526):0.6431909052):0.8489521227,(((t37:0.2137978091,t58:0.7558171267):0.9977718,t17:0.1602738672):0.814617363,t16:0.7358785511):0.5882123332):0.5108074599):0.4033833735,t26:0.4876904353):0.5582309803,((t57:0.01668963907,t59:0.5974151294):0.7298893328,(t47:0.7515056527,t8:0.1309843506):0.8651109848):0.1992067765);";
        } else if (nTaxa == 128) {
            newickTree = "((((((t39:0.7649383547,((t110:0.1021127468,t93:0.7072338203):0.9325676472,t24:0.9118683825):0.2613232969):0.9214477239,((t85:0.9341977697,t70:0.2963797131):0.8458041763,(t46:0.1810089312,t42:0.9887483846):0.5431077231):0.07199716335):0.6049757914,((t62:0.06699732807,(t82:0.9745887131,t53:0.05804512533):0.8596255246):0.2885261178,t8:0.2781210653):0.7050867996):0.8466502437,((t55:0.8496891749,(((((t63:0.4887691024,t9:0.645997202):0.2245028205,t48:0.9381575577):0.4567391849,(t74:0.1225927498,((t128:0.5535472762,t123:0.8024892944):0.2510584455,t118:0.2658069399):0.4912374786):0.03160187067):0.9561984285,((t22:0.7551615268,t109:0.800172599):0.2538118539,(t119:0.1002656904,t1:0.277816873):0.6198285942):0.5669044226):0.698505745,(t16:0.7818064773,t5:0.562224953):0.04990531271):0.6045962132):0.2612644897,(((t37:0.5623362288,t79:0.8031213568):0.3234013342,(t105:0.4027042398,(t89:0.5922396074,t45:0.7514652649):0.7370713451):0.4568462782):0.7229927233,(((t50:0.4704266409,t83:0.1476850766):0.04999543424,(t52:0.1512604994,t10:0.2689684255):0.4731617132):0.9730487335,t96:0.08345450461):0.09669522801):0.3306830805):0.3639469559):0.5413105213,((((t121:0.6780034725,(t78:0.3464668456,(t114:0.4673865528,t54:0.6320368669):0.2723166975):0.4042255948):0.5148424304,(((t68:0.1330636954,(t101:0.3563030539,t47:0.7784714273):0.9826060105):0.7354439551,t29:0.0394046139):0.05661469279,((t61:0.3182016774,t4:0.05380591122):0.5429139391,(t115:0.2914485643,t94:0.8306911949):0.645495665):0.4302359757):0.1157215962):0.9011412384,(((t32:0.3491547701,(t113:0.7938201542,(t14:0.9682615935,t27:0.5906389551):0.9539967934):0.05411513196):0.3683833564,t30:0.2984234286):0.5119977193,(((t125:0.8165114955,t17:0.7278962729):0.1809694394,((t43:0.2378039453,t56:0.7468863048):0.5202735909,t92:0.4784424049):0.7101986918):0.8933412118,((t71:0.5572886774,(t34:0.005902515026,((t58:0.5092027877,t19:0.7256067502):0.5266641052,t64:0.9775454842):0.266974248):0.2581296801):0.9419371316,(t87:0.09344327031,t38:0.7944989395):0.2840698578):0.2243592187):0.2399550427):0.2483833807):0.6911635317,((t51:0.5613704605,(t66:0.1342775486,t127:0.6610968062):0.6256624253):0.841944227,((t72:0.3324858085,(t81:0.7133331513,t7:0.3429733894):0.6819828975):0.5338961948,(((t75:0.5126527266,t36:0.3843118194):0.8487066827,t116:0.8282274823):0.8582237,(t107:0.5139509637,t21:0.846354394):0.7811944324):0.3913210116):0.7906014896):0.8950835387):0.3611667259):0.3451732746,((((t33:0.7840264093,(t23:0.37712236,(t59:0.9087078185,t41:0.9871765384):0.5156125482):0.6490048757):0.1827180537,((((((t65:0.5319485927,(t104:0.6506252484,t12:0.553049742):0.475794496):0.5374843339,(t31:0.4207948151,t111:0.137239418):0.04950388009):0.4507772294,(t84:0.5513098401,t76:0.3503543781):0.7242351691):0.3301062675,t99:0.1935854808):0.2754189558,(t106:0.109507353,t2:0.8231811284):0.9133396847):0.9711034084,((t90:0.9753313179,(t88:0.3354867203,t80:0.9928542636):0.384116922):0.3113528495,(((t44:0.03561024903,(t112:0.8320769672,(t67:0.08186018094,t86:0.6583316473):0.461762442):0.8600317028):0.6091571064,(t13:0.009457860142,t117:0.6073674706):0.5442036549):0.3923571319,(((t18:0.3107851327,((t77:0.341640326,t95:0.3778963678):0.5803708916,t60:0.1119372074):0.7689226302):0.06462831609,((t126:0.2102965228,t100:0.8380059837):0.5803153128,t97:0.9063019322):0.8123584371):0.9084136069,(((t98:0.8817691987,t26:0.2055368342):0.5826433522,(t25:0.7219101344,(t73:0.2072048127,t28:0.5602732815):0.4777206683):0.2312484337):0.5659262224,t3:0.6064614251):0.7041085276):0.1079885087):0.9386998334):0.05993422773):0.4026008353):0.3696106318,((((t6:0.4439277612,t122:0.5753030374):0.924475549,t57:0.5635414056):0.756264915,t102:0.9415287538):0.5020783402,(t120:0.2707241238,(t15:0.1534400461,t108:0.7072247982):0.6599611791):0.8473271159):0.4934334254):0.5243155963,(((t20:0.1073557341,t11:0.1383078767):0.605450385,t69:0.8741544338):0.4238542595,((t40:0.6051922338,t91:0.3390282413):0.03981031221,(((t124:0.5545092474,t49:0.6364520306):0.759151683,t35:0.8348291307):0.1369124381,t103:0.4204189568):0.8130538934):0.3476776269):0.01200653613):0.7374855645);";
        } else {
            throw new IllegalArgumentException("Invalid nTaxa = " + nTaxa);
        }
        System.out.println("Tree is " + newickTree + "\n");
        tree = CodonData.getTree(codonAlignment, newickTree, false);

        System.setProperty("java.only","true");
    }


    // return how many times faster
    private double benchmarking(int iteration) {

        // =============== Standard likelihood ===============
        long[] elapsedTimeMillis = new long[iteration];

        // Get current time
        long start = System.currentTimeMillis();

        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel);

        long standardInit = System.currentTimeMillis()-start;

        double timeStandard = 0;
        for (int i=0; i<iteration; i++) {
            start = System.currentTimeMillis();

            double logP = likelihood.calculateLogP();

            // Get elapsed time in milliseconds
            elapsedTimeMillis[i] = System.currentTimeMillis()-start;
            timeStandard += elapsedTimeMillis[i];

            System.out.println("i = " + i + " logP = " + logP);
        }
        double timeStandardStdev = 0;
        for (int i=0; i<iteration; i++) {
            timeStandardStdev += (elapsedTimeMillis[i]-timeStandard/iteration)*(elapsedTimeMillis[i]-timeStandard/iteration);
        }
        timeStandardStdev = Math.sqrt(timeStandardStdev / iteration);

        // =============== DA likelihood ===============

        int tipCount = tree.getLeafNodeCount();
        int siteCount = codonAlignment.getSiteCount();

        Random generator = new Random(777);
        final int[][] states = generateInternalNodeStatesVertMT(generator, tipCount, siteCount);

        // Get current time
        start = System.currentTimeMillis();

        int internalNodeCount = tree.getInternalNodeCount();
        // vert MT stop codon states = 8, 10, 48, 50
        int[] stopCodons = codonAlignment.getGeneticCode().getStopCodonStates();
        System.out.println("Stop codon states " + Arrays.toString(stopCodons));

        InternalNodeStates internalNodeStates = new InternalNodeStates(internalNodeCount, siteCount);
        // internal nodes
        for (int i=tipCount; i<tree.getNodeCount(); i++)
            internalNodeStates.setNrStates(i, states[i-tipCount]);

        DACodonTreeLikelihood daLikelihood = new DACodonTreeLikelihood();

        daLikelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel,
                "internalNodeStates", internalNodeStates);

        long daInit = System.currentTimeMillis()-start;

        long[] elapsedTimeMillis2 = new long[iteration];
        double timeDA = 0;
        for (int i=0; i<iteration; i++) {
            start = System.currentTimeMillis();

            double logPDA = daLikelihood.calculateLogP();

            // Get elapsed time in milliseconds
            elapsedTimeMillis2[i] = System.currentTimeMillis()-start;
            timeDA += elapsedTimeMillis2[i];

            System.out.println("i = " + i + " DA logP = " + logPDA);
        }
        double timeDAStdev = 0;
        for (int i=0; i<iteration; i++) {
            timeDAStdev += (elapsedTimeMillis2[i]-timeDA/iteration)*(elapsedTimeMillis2[i]-timeDA/iteration);
        }
        timeDAStdev = Math.sqrt(timeDAStdev / iteration);

        // =============== report ===============

        System.out.println("\n=============== Standard likelihood ===============\n");
        System.out.println("Init time " + standardInit + " milliseconds");
        System.out.println(iteration + " iteration(s) " + timeStandard + " milliseconds");
        timeStandard /= iteration;
        System.out.println("Calculation time = " + timeStandard + " +- " + df.format(timeStandardStdev) +
                " milliseconds per iteration in average.");


        System.out.println("\n=============== DA likelihood ===============\n");
        System.out.println("Init time " + daInit + " milliseconds");
        System.out.println(iteration + " iteration(s) " + timeDA + " milliseconds");
        timeDA /= iteration;
        System.out.println("Calculation time = " + timeDA + " +- " + df.format(timeDAStdev) +
                " milliseconds per iteration in average.");

        System.out.println("\n" + df.format(timeStandard/timeDA) + " times faster\n\n");

        return timeStandard/timeDA;
    }

    // vertebrateMitochondrial
    private int[][] generateInternalNodeStatesVertMT(Random generator, int tipCount, int siteCount) {
        System.out.println("Generate internal node states using VertMT : tips = " + tipCount + ", codon = " + siteCount );
        // internal nodes
        int[][] states = new int[tipCount-1][siteCount];

        for (int i=0; i < states.length ; i++) {
//            states[i] = new int[siteCount];
            // 0 - 63
            for (int j=0; j < states[0].length; j++) {
                states[i][j] = (int)(generator.nextDouble() * 64);
                // stop codon states in vertebrateMitochondrial
                while(states[i][j] == 8 || states[i][j] == 10 || states[i][j] == 48 || states[i][j] == 50)
                    states[i][j] = (int)(generator.nextDouble() * 64);
            }
        }
        return states;
    }


    /**
     * Sum 64*64 time is 4577 milliseconds
     * Sum 64 time is 74 milliseconds
     * Take 1 element 4 milliseconds
     *
     * end1 / end3 = 1144 times
     * end2 / end3 = 18 times
     */
    @Test
    public void benchmarkingForLoop(){

        double[][] m = new double[64][64];
        double sum = 0;

        long start = System.currentTimeMillis();
        for (int n=0; n<1000000; n++) {
            for (int i = 0; i < m.length; i++) {
                for (int j = 0; j < m.length; j++) {
                    sum += m[i][j];
                }
            }
        }
        long end1 = System.currentTimeMillis()-start;
        System.out.println("\nSum 64*64 time is " + end1 + " milliseconds");

        start = System.currentTimeMillis();
        for (int n=0; n<1000000; n++) {
            for (int j=0; j<m.length; j++) {
                sum += m[16][j];
            }
        }
        long end2 = System.currentTimeMillis()-start;
        System.out.println("Sum 64 time is " + end2 + " milliseconds");


        start = System.currentTimeMillis();
        for (int n=0; n<1000000; n++) {
            sum += m[16][63];
        }
        long end3 = System.currentTimeMillis()-start;
        System.out.println("Take 1 element " + end3 + " milliseconds");


        System.out.println("\nend1 / end3 = " + (end1 / end3) + " times");
        System.out.println("end2 / end3 = " + (end2 / end3) + " times");
    }

    /**
     * Take reference from 2D array: time is 6 milliseconds
     * System.arraycopy: time is 3562 milliseconds
     * For loop: time is 1806 milliseconds
     */
    @Test
    public void benchmarkingArrays(){
        // 2*16-1
        double[][] nodes = new double[31][10000];
        for (int j = 0; j < nodes[16].length; j++)
            nodes[16][j] = j;
        double[] arr;

        long start = System.currentTimeMillis();
        for (int n=0; n<1000000; n++) {
            arr = nodes[16];
        }
        long end1 = System.currentTimeMillis()-start;
        System.out.println("\nTake reference from 2D array: time is " + end1 + " milliseconds");

        start = System.currentTimeMillis();
        arr = new double[nodes[16].length];
        for (int n=0; n<1000000; n++) {
            System.arraycopy(arr, 0, nodes[16], 0, arr.length);
        }
        long end2 = System.currentTimeMillis()-start;
        System.out.println("System.arraycopy: time is " + end2 + " milliseconds");


        start = System.currentTimeMillis();
        arr = new double[nodes[16].length];
        for (int n=0; n<1000000; n++) {
            for (int j = 0; j < nodes[16].length; j++)
                arr[j] = nodes[16][j];
        }
        long end3 = System.currentTimeMillis()-start;
        System.out.println("For loop: time is " + end3 + " milliseconds");

    }


// ========== in dev ===========

//    private void init6T333(String pi) {
//        Alignment data = CodonTestData.getAlig6T333();
//        // create Codon Alignment
//        codonAlignment = new CodonAlignment();
//        codonAlignment.initByName("data", data, "dataType", "codon",
//                "geneticCode", "vertebrateMitochondrial", "verbose", false);
//
//        if (pi=="F3X4") {
//            String newickTree = "(((Human_Horai: 0.1, Human_Arnason: 0.1): 0.15, (Chimp_Horai: 0.15, " +
//                    "Chimp_Arnason: 0.15): 0.1): 0.25, (Gorilla_Horai: 0.2, Gorilla_Arnason: 0.2): 0.3);";
//
//            CodonFrequencies codonFreq = new CodonFrequencies();
//            codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", false);
//            siteModel = CodonTestData.getSiteModel("0.08000", "15.34858", codonFreq, false);
//            tree = CodonTestData.getTree(codonAlignment, newickTree, false);
//
//            System.setProperty("java.only","true");
//        } else {
//            throw new IllegalArgumentException("Invalid pi " + pi);
//        }
//    }

}
