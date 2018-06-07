"""A parser for reading Gipsy time dependent parameter (TDP) file

Example:
--------

    from where import parsers
    p = parsers.parse_file(parser_name='gipsy_tdp', file_path='final.tdp')
    data = p.as_dict()

Description:
------------

Reads data from files in Gipsy time dependent parameter (TDP) format.


"""

# Where imports
from where.parsers._parser_line import LineParser
from where.lib import plugins


@plugins.register
class GipsyTdpParser(LineParser):
    """A parser for reading Gipsy time dependent parameter (TDP) file
    """

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        TODO: Station name should be separated from parameter name.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file. Following parameters are read:

        ====================  =========================================================================================
         Parameter             Description
        ====================  =========================================================================================
         apriori_value         Nominal value. This field contains the last value used by the model.
         name                  Parameter name.
         sigma                 The sigma associated with the value of the parameter. A negative value indicates it
                               should be used for interpolation by the file reader read_time_variation in
                               $GOA/libsrc/time_variation. If no sigmas are computed by the smapper, a 1.0 will be
                               placed here.
         time_past_j2000       Time given in GPS seconds past J2000.
         value                 Accumulated value of the parameter at time and includes any nominal, or iterative
                               correction. This is the only entry used by the model.
        ====================  =========================================================================================
        """

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----+----8-
        #  564168900.0000   0.00000000000000     -8.944780025556555E-11  5.000E-04  TRPAZSINABSA
        #  564168900.0000   0.00000000000000     -9.807976019950053E-11  5.000E-04  TRPAZCOSABSA
        #  564190800.0000  9.542362619170887E-02  9.542423063851302E-02  0.500      WETZTROPABSA
        #  564190800.0000   0.00000000000000       286.902875647374       388.      STA BIASABSA
        #  564168899.7950   0.00000000000000     -3.338311367575861E-04  -1.00      PB GPS65  ABSA
        #  564168899.8000   0.00000000000000     -3.338311367575861E-04  4.091E-04  PB GPS65  ABSA
        return dict(
            names=("time_past_j2000", "apriori_value", "value", "sigma", "name"),
            delimiter=(15, 23, 23, 11, 20),
            dtype=("f8", "f8", "f8", "f8", "U20"),
            autostrip=True,
        )
