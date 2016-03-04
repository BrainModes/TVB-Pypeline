# coding: utf-8
# Supply some utility functions for the pipeline; various topics!


def mri_convert_bm(in_file, out_file):
    # Use FREESURFERs mri_convert since the Traited interface produces errors!
    from nipype.interfaces.base import CommandLine

    out_orientation = 'RAS'

    cli = CommandLine(command='mri_convert')
    cli.inputs.args = ('--out_orientation ' + out_orientation +
                        ' --input_volume ' + in_file +
                        ' --output_volume ' + out_file)
    cli.run()

    return out_file
