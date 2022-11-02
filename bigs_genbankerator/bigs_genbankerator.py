# This is a script developed by Maddy Butler for batch exporting annotated contigs from private BIGSdb databases 
# It makes use of the BIGSdb RESTful API 
# Currently only works for projects, working on how to adapt it to work for a manual input of isolate ids 

# Importing packages etc. used in script

import click
import os
import datetime 

# for authentication
# used to get session tokens 
# also used to make API calls (rather than requests)
from rauth import OAuth1Service, OAuth1Session

# for sequence manipulations 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation 

# Setting up command line interface using click 
    # mode - if testing, the script will run on a subset of ~10 isolates from the pneumo library
    # if project, the script will run on all the isolates of a specified project, not working yet 

@click.command()
@click.option("--db_url", required=True, help="url of BIGS database to extract annotation files from")
@click.option("--mode", required = True, help="mode to run code, possible values are auth_testing or project")
@click.option("--project", required=False, help="name of project to extract annotation files for, used when mode is project")

def main(db_url, mode, project):
    # use this to make authenticated API calls 
    session_auth = auth(db_url)
    genbank_output(mode, db_url, project, session_auth)

## testing mode deprecated - script is now designed for protected databases only 
## with authentication 

# example database url: "http://rest.pubmlst.org/db/streptococci"
# this is the private streptococci database used by MB when developing this script 

#####____________Defining functions____________#####

##___Authentication function___##
# Returns an authenticated session which is valid for 60 mins 
# Prompts the user to get a verification token from a url 
# This is taken directly from OAuth_testing.py 
# TO DO: ADJUST ALL FUNCTIONS BELOW 

# for accessing streptococci database
# client key and secret provided by KJ 
# consider saving the client key and secret to a file which is read in by the script 

def auth(db_url):
    client_key = "Ckf8Ow12nHNuFuUOPutsEBca"
    client_secret = "Y(L&vm1!U&tgSQZp&XO@hHHQ7IXvOvyz3^pH8sMspq"

    web_url = "https://pubmlst.org/bigsdb?db=streptococci"
    test_web_url = web_url + "&page=authorizeClient"

    request_token_url = db_url + "/oauth/get_request_token"
    access_token_url = db_url + "/oauth/get_access_token"
    session_token_url = db_url + '/oauth/get_session_token'

    # getting the request token 
    service = OAuth1Service(consumer_key=client_key, 
                        consumer_secret=client_secret,
                        request_token_url=request_token_url, 
                        access_token_url = access_token_url, 
                       base_url = test_web_url)
    r = service.get_raw_request_token(params={'oauth_callback':'oob'})

    # using the request token to get the verification url 
    # request token stored in r
    print("Please visit:\n " + web_url + "&page=authorizeClient&oauth_token=" + r.json()["oauth_token"])
    
    # prompting user to add an input
    # the verification token obtained by visiting the above url 
    ver_token = input('Please enter verification code: ')
    
    # getting an access token 
    # using the verification token 
    r2 = service.get_raw_access_token(r.json()["oauth_token"],
                                        r.json()["oauth_token_secret"],
                                        params={'oauth_verifier':ver_token})

    # getting the session token 
    # using the access token stored in r2
    session_request = OAuth1Session(client_key,
                                client_secret, 
                                access_token=r2.json()["oauth_token"], 
                                access_token_secret=r2.json()["oauth_token_secret"])
    r3 = session_request.get(session_token_url)

    # this is an authenticated session 
    session = OAuth1Session(client_key,
                        client_secret, 
                        access_token=r3.json()["oauth_token"], 
                        access_token_secret=r3.json()["oauth_token_secret"])

    print("Authenticated access to database:")
    print(db_url)
    return session 

##___Function to output a new directory full of genbank files for each isolate record within the project___##

def genbank_output(mode, db_url, project, session_auth):

    # create a new directory to save all the genkbank files to 
    date = str(datetime.date.today().strftime("%d%m%Y"))
    
    parent_dir = os.getcwd()
    new_dir = ""
    if mode == "project":
        new_dir = date + "_annotated_contigs_" + str(project)
    elif mode == "auth_testing":
        new_dir = date + "_annotated_contigs_testing"

    path = os.path.join(parent_dir, new_dir)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

    # logging directory generation 
    print("Generated new directory: " + path)

    # call single_isolate_gbk to generate a set of genbank files for a single isolate 

    if mode == "project":
        for isolate_url in get_project_ids(db_url, project, session_auth):
            # isolate_url is generated by get_project_ids 
            # path name is the name of the new directory, tells the function where to write out the genbank files 
            single_isolate_gbk(path, isolate_url, session_auth)
    
    #elif mode == "testing":
    #    for isolate_url in test_isolate_urls(db_url): 
            #isolate_url generated by test_isolate_urls()
    #        single_isolate_gbk(path, isolate_url)

    elif mode == "auth_testing":
        for isolate_url in auth_test_isolate_urls(db_url):
            #isolate_url generated by test_isolate_urls()
            single_isolate_gbk(path, isolate_url, session_auth)

    else:
        print("invalid value given for mode, should be either testing or project")

    # logging after writing out the genbank files  
    print("Generated genbank files for each isolate")

##___Getting a list of isolate id urls associated with the given project___##

def get_project_ids(db_url, project, session_auth):

    # pull the url out as a string from the arguments taken in command line
    # combining the database url and the project name to get a list of isolate id urls 
    
    project_code = "blah"
    if project == "571":
        project_code = "15"
    elif project == "Kenya_Main":
        project_code = "28"
    elif project == "VICE_MaddyDavid":
        project_code = "32"
    elif project == "Non-pneumococcal_species_2021":
        project_code = "33"
    
    if project_code == "blah":
        print("Invalid project name given")
        # need to quit here and re-run with different argument 
    
    url = db_url + "/projects/{}/isolates".format(project_code)

    # returns a list of urls, one for each of the isolates in the specified project 
    # make the api call     
    # include the params so that ALL isolates are returned, not just first 100 
    get_project = session_auth.get(url, params = {"return_all":1}) 
    # check that auth worked
    # status code == 200 when the call has been successful 
    if get_project.status_code == 200:
        pass
    # any other status code indicates something has gone wrong 
    else:
        # try reseting authentication
        # session token may have expired 
        print("getting new session token..")
        session_auth = auth(db_url)
        get_project = session_auth.get(url, params = {"return_all":1})
        if get_project.status_code == 200:
            pass
        else:
            print("Something other than authentication is going wrong")
            print("get_project_ids()")
            # need to quit here 

    # return the list of isolate urls
    print("retrieved isolate url list for project")
    return dict(get_project.json())["isolates"]

    # called within genbank_output(project)

##___Generating a list of isolate id urls for testing the script___##

# command for mode = auth_testing 
def auth_test_isolate_urls(db_url):
    
    # manually defining a list of isolates - one from Iceland, one from Kenya
    # for testing authentication - access to the private streptococcal database 
    id_list = ["322", "6419"]
    
    # populate the list of urls with each isolate from the list
    # manually - not via an API call 
    url_list = []
    for isolate in id_list:
        url_list.append("{}/isolates/".format(db_url) + isolate)
    return url_list


##___Getting a list of contigs urls associated with a given isolate___## 

def get_isolate_contigs(isolate_id_url, session_auth): 

    # generating the url name to access the contigs of the isolate
    url = isolate_id_url + "/contigs"

    # returns a list of urls for contigs 
    #include the params so that ALL contigs are returned, not just first 100 
    get_contigs = session_auth.get(url, params = {"return_all":1}) 
    if get_contigs.status_code == 200:
        pass
    # any other status code indicates something has gone wrong 
    else:
        # try reseting authentication
        # session token may have expired 
        print("getting new session token..")
        session_auth = auth()
        get_contigs = session_auth.get(url, params = {"return_all":1})
        if get_contigs.status_code == 200:
            pass
        else:
            print("Something other than authentication is going wrong")
            print("get_isolate_contigs()")
            print("isolate_id_url")
            # need to quit here 

    return dict(get_contigs.json())["contigs"]

##___Converting a contig record from json to biopython seqrecord format___##

def get_contig_seqrecord(contig_url, session_auth):

    # converting the contig record from json to a dict, easily manipulated to a seqrecord 
    # print(contig_url)

    single_contig_get = session_auth.get(contig_url)
    if single_contig_get.status_code == 200:
        contig_dict = dict(single_contig_get.json())
        return contig_dict
    else:
        session_auth = auth()
        single_contig_get = session_auth.get(contig_url)
        contig_dict = dict(single_contig_get.json())
        return contig_dict

def seqrecord_conversion(contig_dict):

    # generate a list of features, one per annotated locus from the contig record 

    feature_list = []

    # exception written in for contigs without any annotations 
    #     
    if "loci" in contig_dict.keys():
        
        # convert direction of locus from string to numerical
        for locus in contig_dict["loci"]:
            direction = 0
            if locus["direction"] == "forward":
                direction = 1
            elif locus["direction"] == "reverse":
                direction = -1
            else:
                pass

            # generate the seqfeature and populate the empty list 
            # start position was always out by 1, corrected by subtracting 1
            # checked - the output files here now match the BIGS ref export for isolate 1235 

            seq_feature = SeqFeature(location = FeatureLocation(((locus["start"])-1), locus["end"], 
            strand = direction),
            type = "CDS", 
            # id field is not maintained in the genbank exports
            id = locus["locus_name"], 
            # qualifiers gets taken through to the genbank features table
            # work-around for maintaining the name of the loci 
            qualifiers = {"locus":locus["locus_name"]}
            )
            feature_list.append(seq_feature)
    
    else:
        # if no loci annotated, feature_list remains empty 
        pass
    #print(feature_list)

    # manually feed the various fields from the contig record into a SeqRecord object 

    contig_seqrecord = SeqRecord(
        Seq(contig_dict["sequence"]),
        id = str(contig_dict["id"]), 
        features = feature_list, 
        # this argument is required for genbank export, replaces deprecated alphabet attribute to seq
        annotations={"molecule_type": "DNA"}
        )

    # return the contig seqrecord 
    return contig_seqrecord

##___Generating a single genbank format file for all contigs of a given isolate___##

def single_isolate_gbk(dir_path, isolate_url, session_auth):

    isolate_id = isolate_url.split("/")[-1]
    filename = dir_path + "/id_" + isolate_id + "_contigs.gbk" 
    
    # exception - if the file already exists, then don't do it 
    # this should allow me to pick up if the script terminates midway through 
    # e.g. if the connection is lost 
    if os.path.exists(filename):
        print("Contigs already exported for id " + isolate_id)
        return

    seqrecord_list = []

    # Iterates around all the contigs of the given isolate, populates the seqrecord_list

    for contig_url in get_isolate_contigs(isolate_url, session_auth):
        seqrecord_list.append(seqrecord_conversion(get_contig_seqrecord(contig_url, session_auth)))

    # Writes the genkbank file out to the new directory, defined as an argument 
    SeqIO.write(seqrecord_list, filename, "genbank")


#####____________Calling functions____________#####
#print("here")
#genbank_output()

main()