#################################################################
# These global variables have to do with executable locations.
MacLocations = {'build_lmer_table':'/usr/local/RepeatScout/build_lmer_table',
                'RptScout':'/usr/local/RepeatScout/RepeatScout',
                'filter_stage-1':'/usr/local/RepeatScout/filter-stage-1.prl',
                'filter_stage-2':'/usr/local/RepeatScout/filter-stage-2.prl',
                'RAIDER':'./raider',
                'raider_pre':'./raider_pre',
                'bigfoot':'./bigfoot',
                'python':'python3.4',
                'araider':'./araider',
                'phRAIDER': './phRAIDER',
                'pre-phRAIDER':'./pre-phRAIDER',
                'rm_modules': None,
                'RepeatMasker' : 'RepeatMasker',
                'proc_per_node' : 1,
                'basic_arch_type' : None,
                'high_mem_arch' : None, 
                'blast' : 'blastn',
                'blast_modules' : [],
                'time_cmd' : "/usr/bin/time -f \"User: %U\\nM: %M\\nK: %K\"",
                'CompositeDiscover' : './CompositeDiscover',
                'naive' : './naive'}
RedhawkLocations = {'build_lmer_table':'./build_lmer_table',
                    'RptScout':'./RepeatScout',
                    'filter_stage-1':'./filter-stage-1.prl',
                    'filter_stage-2':'./filter-stage-2.prl',
                    'RAIDER':'./raider',
                    'raider_pre':'./raider_pre',
                    'bigfoot':'./bigfoot',
                    'python':'python3.3',
                    'araider':'./araider',
                    'phRAIDER': './phRAIDER',
                    'pre-phRAIDER':'./pre-phRAIDER',
                    'rm_modules' : ['RepeatMasker', 'python-3.3.3'],
                    'RepeatMasker' : 'RepeatMasker',
                    'proc_per_node' : 8,
                    'basic_arch_type' : ["n09","bigmem"],
                    'high_mem_arch' : 'redhawk',
                    'blast' : 'blastn',
                    'blast_modules' : ['blast+'],
                    'time_cmd' : '/usr/bin/time -v',
                    'CompositeDiscover' : './CompositeDiscover',
                    'naive' : './naive'}
OakleyLocations = {'build_lmer_table':'./build_lmer_table',
                   'RptScout':'./RepeatScout',
                   'filter_stage-1':'./filter-stage-1.prl',
                   'filter_stage-2':'./filter-stage-2.prl',
                   'RAIDER':'./raider',
                   'raider_pre':'./raider_pre',
                   'bigfoot':'./bigfoot',
                   'python':'python',
                   'araider':'./araider',
                   'phRAIDER': './phRAIDER',
                   'pre-phRAIDER':'./pre-phRAIDER',
                   'rm_modules' : None,
                   'RepeatMasker' : 'RepeatMasker',
                   'proc_per_node' : 12,
                   'basic_arch_type' : None,
                   'high_mem_arch' : 'oakley',
                   'blast' : 'blastn',
                   'blast_modules' : ['blast'],
                   'time_cmd' : '/usr/bin/time',
                   'CompositeDiscover' : './CompositeDiscover',
                   'naive' : './naive'}



try:
    fp = open("locations.txt")
except: 
    sys.stderr.println("locations.txt file not present.")
    exit(1);

location = fp.readline().strip()
if not (location.lower() in {'redhawk', 'oakley', 'osx'}):
    sys.stderr.println("location.txt: Bad content")
    exit(1)

Locations = None
if location.lower() == "oakley":
    Locations = OakleyLocations;
elif location.lower() == "osx":
    Locations = MacLocations
elif location.lower() == "redhawk":
    Locations = RedhawkLocations
else:
    sys.stderr.println("locations.txt file: bad content")
    exit(1);
