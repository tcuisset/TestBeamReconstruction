import argparse

#custom action to set multiple flags when calling '--all'
class AllAction(argparse.Action):
    def __init__(self, option_strings, dest, vars_to_ignore, **kwargs):
        self.vars_to_ignore = vars_to_ignore
        super(AllAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
         
    def __call__(self, parser, namespace, values, option_string=None):
        for k,v in namespace.__dict__.items():
            if k not in self.vars_to_ignore:
                namespace.__dict__[k] = True
        
class AddArgs():
    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.datatype_choices = ['data', 'sim_proton', 'sim_noproton', 'sim_cmssw']
        self.showertype_choices = ['em', 'had']
    
    def _help_choices(self, s, choices):
        s += ' / '.join(choices)
        return s[:-3] + '.'
    
    def resp_res(self):
        self.parser.add_argument(
            '--analyze_only',
            action='store_true',
            help='Run solely the analysis step.'
        )
        self.parser.add_argument(
            '--plot_only',
            action='store_true',
            help='Run solely the plotting step.'
        )
        
        requiredNamedGroup = self.parser.add_argument_group('required named arguments')
        requiredNamedGroup.add_argument(
            '--datatype',
            type=str,
            choices=self.datatype_choices,
            required=True,
            help=self._help_choices('Choose the datatype to run the analysis on: ', self.datatype_choices)
        )
        requiredNamedGroup.add_argument(
            '--showertype',
            type=str,
            choices=self.showertype_choices,
            required=True,
            help=self._help_choices('Choose the showertype to run the analysis on: ', shower_choices)
        )
        return self.parser.parse_known_args()

    def clusters(self):
        self.parser.add_argument(
            '--hits',
            action='store_true',
            help='Run the cluster analysis on the number of hits per cluster'
        )
        self.parser.add_argument(
            '--energies',
            action='store_true',
            help='Run the cluster analysis on the energy per cluster'
        )
        self.parser.add_argument(
            '--numbers',
            action='store_true',
            help='Run the cluster analysis on the number of clusters'
        )
        self.parser.add_argument(
            '--posx',
            action='store_true',
            help='Run the cluster analysis on the X position of clusters'
        )
        self.parser.add_argument(
            '--posy',
            action='store_true',
            help='Run the cluster analysis on the Y position of clusters'
        )
        self.parser.add_argument(
            '--posx_posy',
            action='store_true',
            help="Run the cluster analysis on the X and Y positions f clusters"
        )
        self.parser.add_argument(
            '--dx',
            action='store_true',
            help="Run the cluster analysis on their X spatial resolution"
        )
        self.parser.add_argument(
            '--dy',
            action='store_true',
            help="Run the cluster analysis on their Y spatial resolution"
        )
        self.parser.add_argument(
            '--dx_dy',
            action='store_true',
            help="Run the cluster analysis on their X and Y spatial resolution (2D plots)"
        )
        self.parser.add_argument(
            '--dx_2D',
            action='store_true',
            help="Run the cluster analysis on their X spatial resolution (2D: resolution vs. layer)"
        )
        self.parser.add_argument(
            '--dy_2D',
            action='store_true',
            help="Run the cluster analysis on their Y spatial resolution (2D: resolution vs. layer)"
        )
        self.parser.add_argument(
            '--use_saved_data',
            action='store_true',
            help='Whether to use the prunned Pandas dataframes stored previously'
        )
        self.parser.add_argument(
            '--chosen_energy',
            type=int,
            default=50,
            help='Beam energy used for all the plots that are layer dependent (positions, resolutions, ...)'
        )

        variables_to_ignore = ['datatype', 'showertype', 'tag']
        self.parser.add_argument(
            '--all',
            action=AllAction,
            vars_to_ignore=variables_to_ignore,
            help='Run the full cluster analysis'
        )
        
        requiredNamedGroup = self.parser.add_argument_group('required named arguments')
        requiredNamedGroup.add_argument(
            '--datatype',
            type=str,
            choices=self.datatype_choices,
            required=True,
            help=self._help_choices('Choose the datatype to run the analysis on: ', self.datatype_choices)
        )
        requiredNamedGroup.add_argument(
            '--showertype',
            type=str,
            choices=self.showertype_choices,
            required=True,
            help=self._help_choices('Choose the showertype to run the analysis on: ', self.showertype_choices)
        )
        requiredNamedGroup.add_argument(
            '--tag',
            type=str,
            required=True,
            help='Identify the tag used for producing the input files.'
        )

        return self.parser.parse_known_args()

    def summary(self):
        self.parser.add_argument(
            '--chosen_energy',
            type=int,
            default=50,
            help='Beam energy used for all the plots that are layer dependent (positions, resolutions, ...)'
        )
        
        requiredNamedGroup = self.parser.add_argument_group('required named arguments')
        requiredNamedGroup.add_argument(
            '--datatype',
            type=str,
            choices=self.datatype_choices,
            required=True,
            help=self._help_choices('Choose the datatype to run the analysis on: ', self.datatype_choices)
        )
        requiredNamedGroup.add_argument(
            '--showertype',
            type=str,
            choices=self.showertype_choices,
            required=True,
            help=self._help_choices('Choose the showertype to run the analysis on: ', self.showertype_choices)
        )
        requiredNamedGroup.add_argument(
            '--var',
            type=str,
            choices=['dx', 'dy'],
            required=True,
            help='Choose the resolution dimension to summarize.'
        )
        
        return self.parser.parse_known_args()

    def layers(self):
        self.parser.add_argument(
            '--densities',
            action='store_true',
            help='Run the layer analysis on the CLUE densities'
        )
        self.parser.add_argument(
            '--distances',
            action='store_true',
            help='Run the layer analysis on the CLUE distances'
        )
        self.parser.add_argument(
            '--densities_distances',
            action='store_true',
            help='Run the layer analysis on the CLUE densities and distances together in the same plot'
        )
        self.parser.add_argument(
            '--densities_2D',
            action='store_true',
            help='Run the layer analysis on the CLUE densities, having all the layer information on the same plot'
        )
        self.parser.add_argument(
            '--distances_2D',
            action='store_true',
            help='Run the layer analysis on the CLUE distances, having all the layer information on the same plot'
        )
        self.parser.add_argument(
            '--hits_fraction',
            action='store_true',
            help='Run the layer analysis on the fraction of clusterized hits'
        )
        self.parser.add_argument(
            '--energy_fraction',
            action='store_true',
            help='Run the layer analysis on the fraction of clusterized energy'
        )
        self.parser.add_argument(
            '--posx_posy',
            action='store_true',
            help="Run the layer analysis on the hits' X and Y positions in the same plot"
        )
        
        variables_to_ignore = ['datatype', 'showertype', 'distances_2D', 'densities_2D', 'posx_posy', 'tag']
        self.parser.add_argument(
            '--all',
            action=AllAction,
            vars_to_ignore=variables_to_ignore,
            help='Run the full layer analysis'
        )
        
        requiredNamedGroup = self.parser.add_argument_group('required named arguments')
        requiredNamedGroup.add_argument(
            '--'+variables_to_ignore[0],
            type=str,
            choices=self.datatype_choices,
            required=True,
            help=self._help_choices('Choose the datatype to run the analysis on: ', self.datatype_choices)
        )
        requiredNamedGroup.add_argument(
            '--'+variables_to_ignore[1],
            type=str,
            choices=self.showertype_choices,
            required=True,
            help=self._help_choices('Choose the showertype to run the analysis on: ', self.showertype_choices)
        )
        requiredNamedGroup.add_argument(
            '--tag',
            type=str,
            required=True,
            help='Identify the tag used for producing the input files.'
        )
        return self.parser.parse_known_args()
