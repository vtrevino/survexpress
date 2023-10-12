<%-- 
    Document   : SurvivaX
    Created on : Nov 10, 2012, 7:11:41 PM
    Author     : trevino
--%>
<%@page import="java.awt.Stroke"%>
<%@ page import="java.util.*, java.io.*"%>
<%@ page import="hibernatePOJOS.*" %> 
<%@ page import="java.sql.ResultSet" %>
<%@ page import="java.sql.SQLException" %>
<%@ page import="java.sql.Statement" %>
<%@ page import="java.text.DateFormat" %>

<%@page contentType="text/html" pageEncoding="UTF-8"%>

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
    "http://www.w3.org/TR/html4/loose.dtd">

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
        <title>SurvExpress - Web resource for Biomarker comparison and validation of Survival gene eXpression data in cancer</title>
        <style type="text/css">
body {
	font-family: verdana,arial,sans-serif;
	background-color: #1B3A59;
	margin-left: 0px;
	margin-top: 0px;
}
.citas {
	font-family: Verdana, Geneva, sans-serif;
	font-size: 10px;
}
table {
	font-family: verdana,arial,sans-serif;
	color: #333333;
	border-width: 1px;
	border-color: #999999;
	border-collapse: collapse;
	background-color: #8ECEFE;
}
a {
	text-align: left;
	font-family: verdana,arial,sans-serif;
	color: #274F87;
}
th {
	background-color: #3A7ABA;
	border-width: 1px;
	padding: 3px;
	border-style: solid;
	border-color: #FFFFFF;
    color: #FFFFFF;
    text-align: left;
    font-size: 14px;
    font-weight:normal;
}
th.function {
	background-color: #FFFFCC;
	border-width: 1px;
	padding: 3px;
	border-style: solid;
	border-color: #a9c6c9;
	text-align: left;
	color: white;
}
td {
	background-color: #FFFFFF;
	border-width: 1px;
	padding: 3px;
	border-style: solid;
	border-color: #a9c6c9;
    font-size: 14px;
	otra-cosa: #b9d6d9;
}
trx:hover {
	border-width: 1px;
	padding: 3px;
	border-style: solid;
	border-color: #a9c6c9;
    background-color:#d4e3e5;
}
.txt {
	font-family: Verdana, Geneva, sans-serif;
	font-size: 12px;
	color: #FFF;
}
th:hover {
	border-width: 1px;
	padding: 3px;
	border-style: solid;
	background-color: #53A4DE;
	color: #FFFFFF;
}
td:hover {
	border-width: 1px;
	padding: 3px;
	border-style: solid;
	background-color: #E7E7E7;
	text-align: left;
}

</style>

<script type="text/javascript" language="javascript">

function replaceString(str, fnd, rep, caseSens) {
	var delim = "/gi";
        if (replaceString.arguments.length > 3) {
            if (caseSens) delim = "/g";
        }
	var regexp = eval("/" + fnd + delim);
	return str.replace(regexp, rep);
}

function replace(st, arr) {
	var i;
	for (i=0; i < arr.length; i++) {
		st = replaceString(st, "%"+(i+1), arr[i], false);
	}
	return st;
}

function replace2(st, arr, caseSens) {
	var i;
	for (i=0; i < arr.length; i++) {
		st = replaceString(st, "%"+(i+1), arr[i], caseSens);
	}
	return st;
}

function changeParamUrlSymbols( url ) {
        url = replaceString(url, "%[^0-9A-Fa-f][^0-9A-Fa-f]","%25", false);
        var symbols = ["'",   "\"", "#",    "&",  "+",   "/",   ":",   ";",   "<",   ">",   "?",    " " ];
        var urls    = ["%27", "%22", "%23", "%26", "%2B", "%2F", "%3A", "%3B", "%3C", "%3E", "%3F", "+" ];
        for (i=0; i < symbols.length; i++) {
                url = url.replace(symbols[i], urls[i]);
        }
        return url;
}

function trim(text) {
    return text.replace(/^\s+|\s+$/g, "");
}

var datasets = new Array();

function AddDataset(key, name, description, authors, options, samples, url, userkey) {
    var a = new Array();
    a[0] = key;
    a[1] = name;
    a[2] = description;
    a[3] = authors;
    a[4] = options.trim().split(",");
    a[5] = samples;
    a[6] = url;
    a[7] = userkey;
    var j;
    for (j=0; j < a[4].length; j++)
        a[4][j] = a[4][j].trim();
    //alert(a[4][0]);
    datasets[datasets.length] = a;
}

var tissueSelected = false;

function changeTissue(tissue) {
    //var radio = document.getElementById("tissue");
    //var tissue = radio.value;
    //alert(tissue);
    var i,j,n=0;
    var s = "<table><th>#</th><th>Database</th><th>Samples</th><th>Clinical data</th><th>Source</th></tr>";
    for (i=0; i < datasets.length; i++) {
        if (datasets[i][4][0] == tissue  ||  tissue == "*") {
            n++;
            s = s + "<tr><td>"+n+"</td><td><input type='radio' id='database"+datasets[i][0]+"' name='database' value='"+datasets[i][0]+"' onclick='changeDataset();'>";
            s = s + "<label for='database"+datasets[i][0]+
                "' title='(Internal ID:"+datasets[i][0]+", UID:"+datasets[i][7]+")\nSummary: "+datasets[i][2]+"\n\nAuthors: "+datasets[i][3]+"'>"+datasets[i][1]+"</label></td><td align=center>" + (datasets[i][5] == 0 ? "" : Math.abs(datasets[i][5]) + (datasets[i][5] < 0 ? "*" : "")) + "</td><td>";
            var datos = "";
            for (j=1; j < datasets[i][4].length; j++) {
                if (j > 1) {
                    datos = datos + ", ";
                }
                datos = datos + datasets[i][4][j];
            }
            a = datasets[i][3].split(/[ ,\t;.]/);
            s = s + "<font size='-" + (datasets[i][4].length > 3 ? datasets[i][4].length-2 : 1) + "'>" + datos + "</font></td>";
            s = s + "<td><a href=\""+datasets[i][6]+"\" target='_blank'>"+a[0]+"</a></td></tr>";
        }
    }
    s = s + "</table>";
    //alert(s);
    document.getElementById("databases").innerHTML = s;
    document.getElementById("send").disabled = true;
    
    return true;
}

    <%
    PrintWriter output = response.getWriter();

    int i;
    File datasets = new File("SurvivaX-datasets.txt");
    
    String ridx = (String) request.getParameter("reindex");
    if ((datasets.exists() ? datasets.lastModified()+86400000*30 : 0) < System.currentTimeMillis() || (ridx != null && "1yestrueok".indexOf(ridx.toLowerCase()) >= 0)) {
        String build = (String) request.getParameter("rebuild");
        boolean rebuild = (build != null && "1yestrueok".indexOf(build.toLowerCase()) >= 0);
        List<DataSet> ld = DataSetDAOHibernate.getPublicDataSetsForBiomarkers();
        if (rebuild) {
            // **** WARNING **** this TAKES A LONG TIME, it is used only for rebuilding index files for each data matrix
            File reindexed = new File("SurvivaX-reindexing.txt");
            FileWriter rfw = new FileWriter(reindexed);
            for (i=0; i < ld.size(); i++) {
                DataSet d = ld.get(i);
                DataMatrix dm = DataMatrixDAOHibernate.getDataMatrixFromDataSetKey(d.getDataSetKey());
                String line = d.getDataSetKey()+"\t"+d.getName()+"\t"+d.getDescription()+"\t"+d.getAuthors()+"\t";
                rfw.write(line+(new Date()).toString());
                rfw.flush();
                if (dm != null) {
                    dm.setRawDataMode(true);
                    dm.buildIndexFile();
                    dm.setRawDataMode(false);
                    dm.buildIndexFile();
                    dm.setRawDataMode(true);
                    rfw.write("\t"+(new Date()).toString()+"\tDONE!\n");
                } else {
                    rfw.write("\t"+(new Date()).toString()+"\tnull matrix\n");
                }
            }
            rfw.close();
        } else {
            FileWriter fw = new FileWriter(datasets);
            for (i=0; i < ld.size(); i++) {
                DataSet d = ld.get(i);
                String options = d.getKeywords();
                options = options.substring(options.lastIndexOf("SurvivaX:")+9);
                if (options.indexOf("/") >= 0) {
                    options = options.substring(0, options.indexOf("/"));
                }
                String[] urls = { d.getAuthorsDataLink(), d.getDataSource(), d.getEbi(), d.getExternalDataLink(), d.getNcbi() };
                String url = null;
                for (int u=0; u < urls.length; u++) {
                    if (urls[u] != null && (url == null || url.length() < urls[u].length())) 
                        url = urls[u];
                }
                String line = d.getDataSetKey()+"^^^"+d.getName()+"^^^"+d.getDescription()+"^^^"+d.getAuthors()+"^^^"+options+"^^^"+url+"^^^"+d.getUserKey();
                line = line.replaceAll("\n|\r|\t", " ").replaceAll("\\^\\^\\^","\t");
                fw.write(line+"\n");
            }
            fw.close();
        }
    }

    Date dsDate = new Date(datasets.lastModified());
    BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(datasets)));
    String line;
    ArrayList<String> tissueList = new ArrayList<String>();
    HashMap<String, String> tissueMap = new HashMap<String, String>();
    int totalSamples = 0;
    while ((line = reader.readLine()) != null) {
        //output.write(line);
        String[] s = line.split("\t");
        String options = s[4];
        int samples = 0;
        int sp = options.toUpperCase().indexOf("SAMPLES");
        if (sp >= 0) {
            String samp = options.substring(sp+7).replaceAll("[^:=#0-9 ]+","");
            try {
                samples = Integer.parseInt(samp.substring(1).replaceAll("[^0-9]+",""));
                if (samp.charAt(0) == '#')
                    samples *= -1;
            } catch (Exception e) {                
            }
            options = options.replaceAll(" ?,?(Samples|SAMPLES|samples).[0-9 ]+", "");
        }
        String[] keywords = options.split(", ?");
        tissueList.add(keywords[0].trim());
        String map = tissueMap.get(keywords[0].trim());
        if (map == null) {
            map = s[1];
        } else {
            map = map+"\n"+s[1];
        }
        tissueMap.put(keywords[0].trim(), map);
        //output.write("AddDataset("+s[0]+",'"+s[1]+"','"+s[2]+"','"+s[3]+"','"+s[4]+"');");
        %>
        AddDataset(<%=s[0].replaceAll("\"|\'|%","")%>,"<%=s[1].replaceAll("\"|\'|%","")%>","<%=s[2].replaceAll("\"|\'|%","")%>","<%=s[3].replaceAll("\"|\'|%","")%>","<%=options.replaceAll("\"|\'|%","")%>",<%=samples%>,"<%=s[5].replaceAll("\"|\'|%","\\\"")%>","<%=s[6].replaceAll("\"|\'|%","\\\"")%>");
        <%
        if (samples > 0)
            totalSamples += samples;
    }
    String tissues[] = tissueList.toArray(new String[0]);
    Arrays.sort(tissues);
    reader.close();    
    
%>
    function setExampleText() {
        var e = document.getElementById("mainform");
        var t = Math.round(Math.random() * (e["tissue"].length-1));
        var tis = e["tissue"][t];
        if (tis != null) {
            tis.click();
        }
        // reference : https://humana-portal.dnadirect.com/pdf/patient-education/humana/oncotypedx.pdf
        //var genes = "MKI67\nESR1\nMMP11\nGRB7\nGSTM1\nACTB\nAURKA\nPGR\nCTSL2\nERBB2\nCD68\nGAPDH\nBIRC5\nBCL2\nBAG1\nRPLP0\nCCNB1\nSCUBE2\nGUSB\nMYBL2\nTFRC".split("\n");
        var genes = "AURKA,BAG1,BCL2,BIRC5,CCNB1,CD68,CTSL2,ERBB2,ESR1,GRB7,GSTM1,MKI67,MMP11,MYBL2,PGR,SCUBE2".split(",");
        var ngenes = Math.round(Math.random() * (genes.length-3))+2;
        var genes2include = "";
        var i;
        for (i=0; i < ngenes; i++) {
            var g = Math.round(Math.random() * (genes.length-1));
            genes2include = genes2include + genes[g] + "\n";
            genes.splice(g, 1);
        }
        alert("A random list of genes will be copied and a random dataset will be selected.\n\nClick the 'SurvExpress Analysis' button when ready.")
        var ds = Math.round(Math.random() * (e["database"].length-1));
        if (e["database"].length == undefined) {
            e["database"].click();
        } else {
            var db = e["database"][ds];
            if (db == null) db = e["database"][0];
            if (db != null) {
                db.click();
            }
        }
        e = document.getElementById("genes");
        e.value = genes2include;
        changeDataset();
        return true;
    }

    function changeDataset() {
        tissueSelected = true;
        enableSurvExpressButton();
        return true;
    }

    function enableSurvExpressButton() {
        var e = document.getElementById("mainform");
        var t = e["tissue"];
        var button = document.getElementById("send");
        var genes = document.getElementById("genes");
        //button.disabled = (((""+genes.value).length == 0) && (t.selectedIndex == undefined));
        button.disabled = (((""+genes.value).length == 0)  ||  (! tissueSelected));
        return true;
    }
    </script>

<script type="text/javascript">
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-37926557-2', 'auto');
  ga('send', 'pageview');

</script>


    </head>
    <body>
        <div style="background: white" align="center">
                            <a href="http://www.itesm.mx" target="_blank"><img src="resources/TecMtyLogo.jpg"></a>
                            <a href="http://www.cmzh.com.mx" target="_blank"><img src="resources/TecSalud.jpg"></a>
                            <a href="http://bioinformatica.mty.itesm.mx" target="_blank"><img src="resources/CatedraBioinformaticaLogo-small.png"></a>
            
        </div>
        <form action="SurvivaXvalidator.jsp" method="post" target="_blank" id="mainform">
        <table>
            <tr>
                <td colspan="2">
                <table>
                    <tr>
                        <td><img src="resources/SurvExpressLogo.png"></td>
                        <td colspan="2" class="citas"><B>Citation:</B><BR>
                            Aguirre-Gamboa R, Gomez-Rueda H, Martínez-Ledesma E, Martínez-Torteya A, Chacolla-Huaringa R, Alberto Rodriguez-Barrientos, José G. Tamez-Peña, Victor Treviño (2013)<I> SurvExpress: An Online Biomarker Validation Tool and Database for Cancer Gene Expression Data Using Survival Analysis. </I><BR>
                            <a target="_blank" href="http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0074250"> PLoS ONE 8(9): e74250. doi:10.1371/journal.pone.0074250</a><BR>
                            <a target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/24066126"> PMID: 24066126 </a> &nbsp;&nbsp;&nbsp;&nbsp;<a target="_blank" href="https://www.dropbox.com/s/pum0zjf0cuu1if1/SurvExpress%20Tutorial%20final-4.pdf?dl=0"> SurvExpress Tutorial </a>
                            <BR>Funding: ITESM grant CAT220, CONACyT grants 83929 and 140601.
                            <BR>Design by Leopoldo Cuéllar + Jesús Abrego</td>
                    </tr>
                </table>
                </td>
            </tr>            
            <tr>
                <th>SurvExpress</th>
                <td>Updated:<br><!-- comment -->
                    - Jul/04/2022 : > 200 Databases re-built<br>
                    - Sep/25/2020 : Link to tutorial<br><!-- comment -->
                    - Oct/xx/2019 : Database access lost & partially recovered<br>
                    - Jan/23/2015 : Faster gene search<br><!-- comment -->
                </td>
            </tr>
            <tr>
                <th>Important Notes</th>
                <td>Given the high demand that occasionally hang the server and the scarce computational and human resources in our lab, this service will be reset each day around 04:00am CST (CST=[UTC/GMT]-6 hours, no summer correction) lasting about one minute.</td>
            </tr>
            <tr>
                <th>(1) Genes:</th>
                <td><select name="geneinputid">
                        <option value="symbol" disabled>Database (Example)</option>
                        <option value="" disabled>------------------</option>
                        <option value="symbol" SELECTED>Symbol (TGFB1)</option>
                        <option value="entrez" >Entrez/GeneID (7040)</option>
                        <option value="Ensembl">Ensembl (ENSG00000153162)</option>
                        <option value="HGNC"   >HGNC (1073)</option>
                        <option value="MIM"    >MIM (112266)</option>
                        <option value="HPRD"   >HPRD (00211)</option>
                        <option value="Vega"   >Vega (OTTHUMG00000014217)</option>
                    </select>
                    <br><textarea name="genes" cols="40" rows="10" id="genes" onkeyup="enableSurvExpressButton();" onchange="enableSurvExpressButton();"></textarea><button type="button" onclick="{ setExampleText(); return false;}">Load Example</button>&nbsp;&nbsp;&nbsp;<a target="_blank" href="http://bioinformatica.mty.itesm.mx/SurvivalToolsTutorials">Tutorial</a></td>
            </tr>
            <tr>
                <th>(2) Tissue:</th>
                <td>
                    <table>
                        <tr>
                    <%
                    String prev = "(nada)";
                    int j=0;                    
                    for (i=0; i < tissues.length; i++) {
                        if ( ! prev.equals(tissues[i].toLowerCase()) ) {
                            String allmaps = tissueMap.get(tissues[i]);
                            if (allmaps == null) allmaps = "";
                            if (j > 0 && (j % 4 == 0)) {
                                %>
                                </tr><tr>
                                <%
                            }
                            //output.write("<td><input id='tissue' name='tissue' type='radio' value='"+tissues[i]+"' onClick='changeTissue()' >"+tissues[i]+"</td>");
                            %>
                            <td>
                                <input id='tissue<%=tissues[i]%>' name='tissue' type='radio' value='<%=tissues[i]%>' onClick='changeTissue("<%=tissues[i]%>")' >
                                <label for='tissue<%=tissues[i]%>' title='<%=allmaps.replaceAll("\"","\'")%>'><%=tissues[i]+" ("+allmaps.split("\n").length+")"%></label>
                            </td>
                            <%
                            j++;
                            prev = tissues[i].toLowerCase();
                        }
                    }
                    if (j > 0 && (j % 4 == 0)) {
                        %>
                        </tr><tr>
                        <%
                    }
                    %>
                            <td>
                                <input id='tissueALLTISSUES' name='tissue' type='radio' value='all' onClick='changeTissue("*")' >
                                <label for='tissueALLTISSUES' title='List all datasets.'>All (<%=tissues.length%>)</label>
                            </td>
                        </tr>
                        <tr>
                            <th colspan="4">Notes:</th>
                        </tr>
                        <tr>
                            <td colspan="4">Tissue or preferred database not listed ? (or found an error) <br>Please, share your data with us by e-mailing corresponding author.</td>
                        </tr>
                        <tr>
                            <td colspan="4">For TCGA datasets, please see <a href=" http://cancergenome.nih.gov/abouttcga/policies/publicationguidelines" target="_blank">TCGA Publication Guidelines</a></td>
                        </tr>
                        <tr>
                            <td colspan="4">Please cite SurvExpress and datasets authors properly.</td>
                        </tr>
                    </table>
                </td>
            </tr>
            <tr>
                <th>(3) Database:</th>
                <td>
                    <div id="databases">Select a tissue</div>
                </td>
            </tr>
            <tr>
                <th>(4) Options:</th>
                <td>
                    <table>
                        <tr>
                            <th>(a) Duplicated genes:</th>
                            <td> This applies when a gene is associated to many rows of the dataset.<br>For example, when a gene has several probe sets (duplicates or alternatives).<hr>
                                <input id="dupmean" name="duplicates" type="radio" value="mean" checked>
                                <label for="dupmean">Average : All probe sets/records will be averaged per sample.</label><br>
                                <input id="dupmax" name="duplicates" type="radio" value="max">
                                <label for="dupmax">Maximum average : The "most expressed" row will be used.</label><br>
                                <input id="dupsd" name="duplicates" type="radio" value="var">
                                <label for="dupsd">Maximum variance : The "most dispersed" row will be used.</label><br>
                                <input id="dupnone" name="duplicates" type="radio" value="all">
                                <label for="dupnone">Show all : All rows will be presented in the analysis. Use this if you have no clue.</label><br>
                            </td>
                        </tr>
                        <tr>
                            <th>(b) Data:</th>
                            <td>
                                <input id="sourceoriginal" name="datasource" type="radio" value="raw" checked>
                                <label for="sourceoriginal">Original (Quantile-Normalized)</label><br>
                                <input id="sourceunif" name="datasource" type="radio" value="uniformized">
                                <label for="sourceunif">Uniformized</label>
                            </td>
                            </tr>
                    </table>
                </td>
            </tr>
            <tr>
                <th>(5) Send:</th>
                <td><input id="send" name="send" type="submit" value="SurvExpress Analysis" disabled>
                </td>
            </tr>
        </table>
        </form>
    <SCRIPT type="text/javascript" language="javascript">
        //changeTissue("");
    </SCRIPT>
    </body>
</html>
