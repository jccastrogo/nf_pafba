#!/usr/bin/env julia
###########################################################################
# @name: create_metabolic_DB.jl                                           #
# @author: Juan C. Castro <jccastrog at gatech dot edu>                   #
# @update: 20-Apr-2018                                                    #
# @version: 1.0.6                                                         #
# @license: GNU General Public License v3.0.                              #
# please type "./create_metabolic_DB.jl -h" for usage help                #
###########################################################################

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
using MySQL;
using ConfParser;
# 1.2 Define funcions ====================================================#
"""
    parse_commandline()

Get the input arguments from ARGV Open the database and insert, modify,
or remove reactions or compounds.

# Examples
```julia-repl
julia> args = parse_commandline();
1
```
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--reactions_file", "-r"
            help = "A table  with metabolic reactions"
			required = true
        "--compounds_file", "-c"
            help = "A table with metabolic compounds"
			required = true
		"--database_conf", "-d"
			help = "A configuraton file for the database"
			required = true
        "--index"
            help = "If present reactions and compounds will be added to the database"
            action = :store_true
    end

    return parse_args(s)
end
"""
    manipulate_db(configFile, operation, [values])

Open the database and insert, modify, or remove reactions,
compounds, or genes. Operations can be any of the following


# Examples
```julia-repl
julia> manipulate_db("configFile.txt", "insert_reaction", values);
1
```
"""
function manipulate_db(configFile, operation, values)
    conf = ConfParse(configFile);
    parse_conf!(conf);
    user = retrieve(conf, "database", "USER");
    password = retrieve(conf, "database", "PASSWORD");
    port = retrieve(conf, "database", "PORT");
    host = retrieve(conf, "database", "HOST");
    socket = retrieve(conf, "database", "SOCKET");
    db=retrieve(conf, "database", "DB");
    host = String(host);
    user = String(user);
    password = String(password);
    db = String(db);
    port = parse(Int64, port)
    conn = MySQL.connect(host, user, password; db=db, port=port, opts = Dict());
    if operation=="insert_reaction"
        reaction_insert(conn, values)
    end
    if operation=="insert_gene"
        gene_insert(conn, values)
    end
    MySQL.disconnect(conn::MySQL.Connection)
end
"""
    reaction_insert(values)

Execute MySQL insert into the database for a set of values.
Values is an array than includes all of the following fields:
id
direction
compartment
name
deltaG
reference
equation
equation_metnames
bigg_id
kegg_id
kegg_pathway
metacyc_pathway

If a value is missing can be set to "NULL".

# Examples
```julia-repl
julia> reaction_insert(["id", "dir", "comp", "name", "dG", "ref", "eq", "eq_met", "bigg_id", "kegg_id", "kegg_pw", "metacyc_pw"]);
1
```
"""
function reaction_insert(conn, values)
    sqlInsert = "INSERT INTO reactions (id, direction, compartment, name, deltaG, reference, equation, equation_metnames, bigg_id, kegg_id, kegg_pathway, metacyc_pathway)";
    sqlValues = values;
    sql = "$sqlInsert $sqlValues";
    try
        MySQL.execute!(conn, sql);
    catch
        write(STDOUT, "WARNING! Could not insert reaction into database \n")
        write(STDOUT, "$sqlValues \n")
    end

end
"""
    gene_insert(values)

Execute MySQL insert into the database for a set of values.
Values is an array than includes all of the following fields:
id
gene

# Examples
```julia-repl
julia> gene_insert(["id", "gene"]);
1
```
"""
function gene_insert(conn, values)
    sqlInsert = "INSERT INTO genes (id, gene)";
    sqlValues = values;
    sql = "$sqlInsert $sqlValues"
    try
        MySQL.execute!(conn, sql);
    catch
        #write(STDOUT, "WARNING! Could not insert gene into database \n")
        #write(STDOUT, "$sqlValues \n")
        write(STDOUT, "$sql \n")
    end

end
"""
    parse_reaction_insert(reactionFile, configFile)

Open the reactions file parse the reaction fields
and add missing reactions to the database.

# Examples
```julia-repl
julia> parse_reaction_file("reactions.txt", "configFile.txt");
1
```
"""
function parse_reaction_insert(reactionFile, configFile)
    rFile = open(reactionFile);
    for line in eachline(rFile)
        line = rstrip(line);
        fields = split(line,"\t");
        id = split(fields[1], "_")[1];
        direction = fields[2];
        compartment = fields[3];
        gpr = fields[4];
        name = fields[5];
        enzyme = fields[6];
        deltag = fields[7];
        reference = "NULL";
        equation = fields[9];
        equationMetnames = fields[10];
        biggID = "NULL";
        keggID = "NULL";
        keggPathway = "NULL";
        metacycPathway = "NULL";
        #Modify values
        name = replace(name, "'", "");
        if fields[7]==""
                deltag = "NULL";
        end
        try
            if fields[8]==""
                reference = "NULL";
            else
                reference = fields[8];
            end
        catch
            reference = "NULL";
        end
        equationMetnames = replace(equationMetnames, "'", "");
        try
            biggID = fields[11]
        catch
            biggID = "NULL";
        end
        try
            keggID = fields[12]
        catch
            keggID = "NULL";
        end
        try
            keggPathways = fields[13]
        catch
                keggPathways = "NULL";
        end
        try
            metacycPathway = fields[14]
        catch
            metacycPathway = "NULL";
        end
        sqlValues = "VALUES ('$id', '$direction', '$compartment', '$name', $deltag, '$reference', '$equation', '$equationMetnames', '$biggID', '$keggID', '$keggPathway', '$metacycPathway');";
        manipulate_db(configFile, "insert_reaction", sqlValues)
    end
    close(rFile)
end
"""
    parse_gene_insert(reactionFile, configFile)

Open the reactions file parse the reaction fields
and add missing reactions to the database.

# Examples
```julia-repl
julia> parse_reaction_file("reactions.txt", "configFile.txt");
1
```
"""
function parse_gene_insert(reactionFile, configFile)
    rFile = open(reactionFile);
    for line in eachline(rFile)
        line = rstrip(line);
        fields = split(line,"\t");
        id = split(fields[1], "_")[1];
        gprArray = fields[4];
        gprArray = replace(gprArray, "(", "")
        gprArray = replace(gprArray, ")", "")
        gprArray = split(gprArray, " and ")
        for gpr in gprArray
            geneArray = split(gpr, "or")
            for gene in geneArray
                rxn = id
                sqlValues = "VALUES ('$rxn', '$gene');";
                manipulate_db(configFile, "insert_gene", sqlValues)
            end
        end
    end
    close(rFile)
end
###========================= 2.0 Main function =========================###
"""
    main()

Execute the functions in order to parse both reactions and
compounds files.

# Examples
```julia-repl
julia> main();
1
```
"""
function main()
    parsedArgs = parse_commandline();
    write(STDOUT, "Parsing reactions file\n");
    reactionFile = parsedArgs["reactions_file"];
    compounds_file = parsedArgs["compounds_file"]
    databaseConf = parsedArgs["database_conf"];
	#parse_reaction_insert(reactionFile, databaseConf);
    parse_gene_insert(reactionFile, databaseConf);
end

###===================== 3.0 Execute main function =====================###
main()
###=====================================================================###

