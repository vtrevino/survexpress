<%-- 
    Document   : SurvivaXvalidator
    Created on : Nov 11, 2012, 12:37:53 AM
    Author     : trevino
--%>
<%@page import="java.text.SimpleDateFormat"%>
<%@page import="java.text.DateFormat"%>
<%@page import="biomatec.SessionBean1" %>
<%@page import="java.util.*, java.io.*" %>
<%@page import="hibernatePOJOS.*" %> s
<%@page import="java.sql.ResultSet" %>
<%@page import="java.sql.SQLException" %>
<%@page import="java.sql.Statement" %>
<%@page import="tools.*" %>

<%@page contentType="text/html" pageEncoding="UTF-8"%>
<!DOCTYPE html>

<%!

static int requestCount = 0;
static DateFormat dateFormat = new SimpleDateFormat("[yyyy/MM/dd HH:mm:ss] ");
static FileWriter log = null;
static int logCount = 0;

static String[] infoLines = null;
static int[] idxColumns = new int[] { 1, 2, 4, 5, 7 };
static int[] columnToIndex = new int[] { -1, 0, 1, -1, 2, 3, -1, 4 };
static Vector<HashMap<String, Integer>> hmKeys = new Vector<HashMap<String, Integer>>(5);
static boolean cached = loadGeneInfoFile();

static boolean loadGeneInfoFile() {
    
    Vector<String> lr = new Vector<String>(50000);
    String line;
    String fields[];
    String keys[];
    int i,j,k;
    for (i=0; i < idxColumns.length; i++) {
        hmKeys.add(new HashMap<String, Integer>(50000));
    }
    try {
        BufferedReader bri = new BufferedReader(new FileReader("/biomatecdatasets/applications/Homo_sapiens.gene_info"));
        while ((line = bri.readLine()) != null) {
            lr.add(line);
        }
        bri.close();
        infoLines = lr.toArray(new String[0]);
        lr = null;
        for (i=0; i < infoLines.length; i++) {
            Integer iLine = new Integer(i);
            fields = infoLines[i].split("\t");
            int min = Math.min(idxColumns.length, fields.length);
            for (j=0; j < min; j++) {
                keys = fields[j].split("\\|");
                HashMap<String, Integer> hm = hmKeys.get(j);
                for (k=0; k < keys.length; k++) {
                    hm.put(keys[k].toUpperCase(), iLine);
                }
            }
        }
    } catch (Exception e) {
        System.out.println("Error reading Homo_sapiens.gene_info.");
        Log("Error reading Homo_sapiens.gene_info.", -1);
        System.out.println(e.getMessage());
        e.printStackTrace();
        return false;
    }
    
    return true;
}


static void Log(String s, int req) {
    Calendar cal = Calendar.getInstance();
    try {
        if (log == null)
            log = new FileWriter("/biomatecdatasets/datasets/SurvivaXvalidator-log.txt", true);
        log.write(dateFormat.format(cal.getTime())+"["+req+"] ");
        log.write(s);
        log.write("\n");
        if (++logCount % 10 == 0) 
            log.flush();
    } catch (Exception e) {
        System.out.println(e.getMessage());
    }
}

static String[] getPossibleGenes(String gene, int field, int output, int req) {
 
    if (!cached) {
        Log("Loading Gene Info.", req);
        cached = loadGeneInfoFile();
        Log("Gene Info. Loaded.", req);
    }
    if (cached && gene.indexOf("*") < 0) {
        Integer l = hmKeys.get(columnToIndex[field-1]).get(gene.toUpperCase());
        if (l != null) {
            String fields[] = infoLines[l].split("\t");
            return new String[] { fields[output-1] };
        }
    }
    
    try {
        Log("find_gene_symbol.sh for ["+gene+"]", req);
        java.lang.Process p = Runtime.getRuntime().exec(new String[] {
            "bash",
            "find_gene_symbol.sh",
            gene,
            ""+field,
            ""+output
        });
        String line;
        String lines="";
        Thread.sleep(100);
        BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
        while ((line = bri.readLine()) != null) {
            if (lines.length() > 0) lines += "\n";
            lines += line; 
        }
        bri.close();
        bri = new BufferedReader(new InputStreamReader(p.getErrorStream()));
        while ((line = bri.readLine()) != null) {
            if (lines.length() > 0) lines += "\n";
            lines += line; 
        }
        bri.close();
        p.getOutputStream().close();
        p.waitFor();
        if (lines.length() == 0) return null;
        return lines.split("\n");
    } catch (Exception e) {
        Log("Error processing bash command:"+e.getMessage(), req);
        return null;
    }

}
%>
<%

requestCount++;
int requestId = requestCount;

int nParams = 0;
Enumeration<String> pnames = request.getParameterNames();
while (pnames.hasMoreElements()) {
    String p = pnames.nextElement();
    Log("Parameter: ["+p+"] = ["+request.getParameter(p).replaceAll("\n|\r|\f",",").replaceAll(",,",",")+"]", requestId);
    nParams++;
}
if (nParams < 4) {
    %>
    <html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
        <title>Biomarker Viewer Error</title>
    </head>
    <body>
     ERROR: NO PARAMETERS
    </body>
    </html>
    <%
    return;
}

String geneInput = request.getParameter("genes");
long dsk = Long.parseLong(request.getParameter("database")); // 
String duplicates = request.getParameter("duplicates"); // mean, max, all
String datasource = request.getParameter("datasource"); // raw, uniformized
String geneID = request.getParameter("geneinputid");

//Log("Getting session");
SessionBean1 sb1 = (SessionBean1) request.getSession(true).getAttribute("SessionBean1");
if (sb1 == null) {
    sb1 = new SessionBean1();
    request.getSession().setAttribute("SessionBean1", sb1);
}
sb1.removeAllSelections();
//Log("Getting user");
User u = UserDAOHibernate.getUserFromKey(new Long(8)); // 8 is the guest 
sb1.setIsGuest(false);
sb1.setLoggedIn(true);
sb1.setLoggedUserId(u.getUserId());
sb1.setLoggedUser(u);
sb1.setIsUserAdmin(false);
//Log("Loading sripts");
sb1.loadScripts();
//Log("Getting dataset");
DataSet ds = DataSetDAOHibernate.getDataSetFromKey(dsk);
//Log("Getting data matrix");
DataMatrix matrix = DataMatrixDAOHibernate.getDataMatrixFromDataSetKey(dsk);
sb1.setEditingDataSet(ds);
sb1.setEditingDataMatrix(matrix);
sb1.message = "Biomarker: Cox Survival Analysis\nSurvExpress - Web resource for Biomarker comparison and validation of Survival gene eXpression data. <BR> Dataset:"+ds.getName()+" (dup="+duplicates+", data="+datasource+")@@GNF@@\\n";

String[] genes = geneInput.split("\n|\r|\f| +|,|;|\\|");

matrix.setRawDataMode(datasource.equals("raw"));
matrix.setTempBufferLinesToSmall(); // to speed up reading

int i,j,k;
String colTypes = matrix.getColumnsType();
String possibleGenes[];
Log("Symbol Matrix File: "+matrix.getIndexPhysicalFileFromUploaded(matrix.isRawDataMode()), requestId);
String genesNotFound = "";
for (i=0; i < genes.length; i++) {
    String gene = genes[i].trim().toUpperCase();
    String inputGene = gene;
    String searchAs = gene;
    if (gene.length() > 0) {
        int field = 0;
        if (geneID.equals("symbol")) {
        } else if (geneID.equals("entrez")) {
            field = 2;
        } else if (geneID.equals("Ensembl")) {
            field = 6;
            if (gene.indexOf("ENSEMBL:") < 0) gene = "Ensembl:"+gene;
        } else if (geneID.equals("HGNC")) {
            field = 6;
            if (gene.indexOf("HGNC:") < 0) gene = "HGNC:"+gene;
        } else if (geneID.equals("MIM")) {
            field = 6;
            if (gene.indexOf("MIM:") < 0) gene = "MIM:"+gene;
        } else if (geneID.equals("HPRD")) {
            field = 6;
            if (gene.indexOf("HPRD:") < 0) gene = "HPRD:"+gene;
        } else if (geneID.equals("Vega")) {
            field = 6;
            if (gene.indexOf("VEGA:") < 0) gene = "VEGA:"+gene;
        }
        if (field != 0) {            
            searchAs = gene;
            possibleGenes = getPossibleGenes(gene, field, 3, requestId);
            if (possibleGenes == null || possibleGenes.length == 0) {
                //gene = "";
                searchAs += "("+geneID+")";
            } else {
                if (possibleGenes.length == 1) {
                    searchAs += "("+geneID+")";
                    gene = possibleGenes[0].trim().toUpperCase();
                } else {
                    searchAs += "( could be ";
                    for (int pg=0; pg < possibleGenes.length; pg++) {
                        searchAs += (pg > 0 ? "," : "") + possibleGenes[pg];
                    }
                    searchAs += ")";
                }
            }
        }
        boolean exact = ! gene.endsWith("*");
        gene = gene.replaceAll("[^A-Z0-9-]","");
        if (gene.length() > 0) {
            if (!gene.equals(inputGene)) 
                searchAs += (searchAs.length() > 0 ? ", " : "") + gene;
            Vector<Integer> vlines = matrix.searchSymbolInIndexFile(gene, false, exact); // true
            Log("Search 1 - Lines "+(exact ? "exact" : "inexact")+" for symbol ["+inputGene+"/"+gene+"] = "+(vlines == null ? "null" : ""+vlines.size()), requestId);
            if ((vlines == null || vlines.size() == 0) && (field==0  ||  field==2)) {
                // NOT FOUND, look for synonyms and return the official symbol
                possibleGenes = getPossibleGenes(gene, 5, 3, requestId);
                if (possibleGenes == null || possibleGenes.length == 0) {
                    //gene = "";
                } else {
                    if (possibleGenes.length == 1) {
                        String newGene = possibleGenes[0].trim().toUpperCase();
                        if (! gene.equals(newGene)  &&  newGene.length() > 1) {
                            //gene = newGene;
                            vlines = matrix.searchSymbolInIndexFile(newGene, false, exact); // search official symbol : true
                            searchAs += (searchAs.length() > 0 ? ", " : "") + newGene;
                            Log("Search 2 - Lines "+(exact ? "exact" : "inexact")+" for symbol ["+inputGene+"/"+newGene+"] = "+(vlines == null ? "null" : ""+vlines.size()), requestId);
                        }
                        String foundAs = null;
                        String sougth = "";
                        if ((vlines == null || vlines.size() == 0)) {
                            // Search data for all synonym symbols
                            possibleGenes = getPossibleGenes(gene, 5, 5, requestId);
                            Vector<Integer> allLines = new Vector<Integer>();
                            for (int ig=0; ig < possibleGenes.length; ig++) {
                                String[] allGenes = possibleGenes[ig].split("\n|\\|");
                                for (int jg=0; jg < allGenes.length; jg++) {
                                    String thisGene = allGenes[jg].trim().toUpperCase();
                                    if (! gene.equals(thisGene) &&  thisGene.length() > 1) {
                                        sougth = (sougth.length() == 0 ? "" : sougth + ", ") + thisGene;
                                        Vector<Integer> rows = matrix.searchSymbolInIndexFile(thisGene, false, exact); // : true
                                        searchAs += (searchAs.length() > 0 ? ", " : "") + thisGene;
                                        Log("Search 3 - Lines "+(exact ? "exact" : "inexact")+" for symbol ["+inputGene+"/"+thisGene+"] = "+(rows == null ? "null" : ""+rows.size()), requestId);
                                        if (rows != null && rows.size() > 0) {
                                            allLines.addAll(rows);
                                            if (foundAs != null) {
                                                foundAs = foundAs + "|" + thisGene;
                                            } else {
                                                foundAs = thisGene;                                            
                                            }
                                        }
                                    }
                                }
                            }
                            if (foundAs == null) {
                                //if (sougth.length() > 0) tag = tag + " Neihter found as "+sougth;
                            } else if (foundAs.indexOf("|") > 0) {
                                //tag = tag + " Warning: found for aliases-symbols "+foundAs.replaceAll("\\|",", ")+".";
                                searchAs += " (Warning: found for "+foundAs.replaceAll("\\|",", ")+")";
                                vlines = allLines;
                            } else {
                                //tag = tag + " Found for alias-symbol "+foundAs+".";
                                searchAs += " (found for "+foundAs+")";
                                vlines = allLines;
                            }
                        }  
                   } else {                 
                        searchAs += " (Ambiguity, possibles:";
                        for (int pg=0; pg < possibleGenes.length; pg++) {
                            searchAs += (pg > 0 ? "," : "") + possibleGenes[pg];
                        }
                        searchAs += ")";
                    }
                }
            }
            if (vlines == null || vlines.size() == 0) {
                genesNotFound = genesNotFound + (genesNotFound.length() > 0 ? ", " : "") + inputGene;
            }
            sb1.message += "\\n" + inputGene + ":" + (vlines == null ? 0 : vlines.size()) + " row(s)." +
                           (!inputGene.equals(gene) ? " Gene:" + gene + "." : "") +
                           (!inputGene.equals(searchAs) ? " Sought as " + searchAs + "." : "");
            if (vlines != null && vlines.size() > 0) {
                if (duplicates.equals("all")) {
                    //for (j=0; j < vlines.size(); j++) {
                    //    sb1.addDsSelectedRows(vlines.get(j)-matrix.getDataHeaderLines());
                    //}
                    Vector<String> sv = matrix.getLinesAtPointer(Utilities.Integer2int(vlines.toArray(new Integer[0])));
                    for (j=0; j < vlines.size(); j++) {
                        sb1.addDsSelectedText(sv.get(j));
                    }
                } else if (vlines.size() == 1) {
                    //sb1.addDsSelectedText(matrix.getLine(vlines.get(0)));
                    sb1.addDsSelectedText((matrix.getLinesAtPointer(new int[] { vlines.get(0).intValue() })).get(0));
                } else if (duplicates.equals("mean")) {
                    Vector<String> sv = matrix.getLinesAtPointer(Utilities.Integer2int(vlines.toArray(new Integer[0])));
                    String geneCols[][] = new String[vlines.size()][];
                    for (j=0; j < vlines.size(); j++) {
                        //geneCols[j] = matrix.getLine(vlines.get(j)).split("\t");
                        geneCols[j] = sv.get(j).split("\t");
                    }
                    for (k=0; k < colTypes.length(); k++) {
                        if (colTypes.charAt(k) == 'D') {
                            double sum = 0;
                            int n = 0;                    
                            for (j=0; j < vlines.size(); j++) {
                                try {
                                    double v = Double.parseDouble(geneCols[j][k]);
                                    sum += v;
                                    n++;
                                } catch (Exception e) {

                                }
                            }
                            if (n > 0) 
                                geneCols[0][k] = ""+Math.round(sum*100000/n)/((double) 100000);
                        } else {
                            /* For survexpress we will not include this to avoid problems in the figures
                            String values = (geneCols[0].length > k ? geneCols[0][k] : "");
                            String raws = values + "$$$";
                            int others = 0;
                            for (j=1; j < vlines.size(); j++) {
                                if (geneCols[j].length > k && raws.indexOf(geneCols[j][k]+"$$$") < 0) {
                                    values = values + ", " + geneCols[j][k];
                                    raws = raws + geneCols[j][k] + "$$$";
                                    others++;
                                }
                            }
                            if (geneCols[0].length > k) {
                                //geneCols[0][k] = values;
                                if (others > 0) {
                                    geneCols[0][k] = geneCols[0][k] + " ["+(others+1)+"]";
                                }
                            }
                            */
                        }
                    }
                    String line = geneCols[0][0];
                    int mxLen = Math.min(geneCols[0].length, colTypes.length());
                    for (k=1; k < mxLen; k++) {
                            line = line + "\t" + geneCols[0][k];
                    }
                    sb1.addDsSelectedText(line);
                } else if (duplicates.equals("var") || duplicates.equals("max")) {
                    Vector<String> sv = matrix.getLinesAtPointer(Utilities.Integer2int(vlines.toArray(new Integer[0])));
                    String geneCols[][] = new String[vlines.size()][];
                    for (j=0; j < vlines.size(); j++) {
                        //geneCols[j] = matrix.getLine(vlines.get(j)).split("\t");
                        geneCols[j] = sv.get(j).split("\t");
                    }
                    StatCalc sc[] = new StatCalc[vlines.size()];
                    for (j=0; j < vlines.size(); j++) {
                        sc[j] = new StatCalc();
                        for (k=0; k < colTypes.length(); k++) {
                            if (colTypes.charAt(k) == 'D') {
                                try {
                                    sc[j].enter(Double.parseDouble(geneCols[j][k]));
                                } catch (Exception e){
                                }
                            }
                        }
                    }
                    int max = 0;
                    for (j=1; j < vlines.size(); j++) {
                        double vmax = (duplicates.equals("var") ? sc[max].getVariance() : sc[max].getMean());
                        double vj   = (duplicates.equals("var") ? sc[j  ].getVariance() : sc[j  ].getMean());
                        if (vj > vmax) {
                            max = j;
                        }
                    }
                    String gc[] = geneCols[max];
                    String line = gc[0];
                    int gclen = gc.length;
                    for (k=1; k < colTypes.length(); k++) {
                        line = line + "\t" + (gclen > k ? gc[k] : "");
                    }
                    sb1.addDsSelectedText(line);
                }
            }
        }
    }
}
Log("Forwarding to BiomarkerViewer.jsp", requestId);
if (log != null) {
    log.flush();
}
log = null;
sb1.message = sb1.message.replaceAll("@@GNF@@", "\\\\n" + (genesNotFound.length() > 0 ? genesNotFound.split(",").length + " Genes not found: "+genesNotFound : ""));
sb1.message += "\\n\\nNext steps: Select censored option, select stratification (optional), click Go.";
%>
<jsp:forward page="BiomarkerViewer.jsp" />

