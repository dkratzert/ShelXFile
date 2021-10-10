# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
from __future__ import print_function

import os
import re
import subprocess
import sys
from shutil import which, disk_usage, copyfile

from shelxfile.misc.misc import remove_file, sep_line, find_line
from shelxfile.shelx.cards import ACTA


def get_xl_version_string(exe: str) -> str:
    """
    Extracts the version string from a SHELXL executable.
    This is fast and needs no hashes etc.
    :type exe: str
    :param exe: path to SHELXL executable
    """
    try:
        with open(exe, 'rb') as f:
            binary = f.read()
            position = binary.find(b'Version 201')
            if position > 0:
                f.seek(position + 8, 0)  # seek to version string
                version = f.read(6)  # read version string
                return version.decode('ascii')
            else:
                return ''
    except IOError:
        print("Could not determine SHELXL version. Refinement might fail to run.")
        return ''


class ShelxlRefine():
    """
    A class to do a shelxl refinement. It is only for shelxl 2017 and above!
    The resfilename should be without ending.
    """

    def __init__(self, shx, resfile_name: str, shelxpath: str = None):
        self.shx = shx
        self.shelxpath = shelxpath
        self.resfile_name, ext = os.path.splitext(resfile_name)
        self._shelx_command = self.find_shelxl_exe(shelxpath)
        self.backup_file = os.path.abspath(str(self.resfile_name + '.shx-bak'))
        self._acta_card = ''  # stores the ACTA values if acta is removed before refinement

        if not self._shelx_command:
            print('\nSHELXL executable not found in system path.\n')
            print('You can download SHELXL at http://shelx.uni-goettingen.de\n')

    def find_shelxl_exe(self, shelxpath=None) -> str:
        """
        returns the appropriate shelxl executable
        """
        names = ['shelxl', 'xl']
        download = 'You can download SHELXL at http://shelx.uni-goettingen.de'
        exe = ''
        if shelxpath:
            return shelxpath
        for name in names:
            exe = which(name)
            if not exe:
                continue
            version = get_xl_version_string(exe)
            if not version:
                print('Your SHELXL version', exe, 'is too old for this Program')
                print('Please use SHELXL 2017 or above!')
                print(download)
            version = version.split('/')
            if int(version[0]) < 2017:
                print('Your SHELXL version is too old. Please use SHELXL 2017 or above!')
                print(download)
            else:
                return exe
        return exe

    def get_b_array(self):
        """
        Approximates the B array size to ensure refinement.
        """
        number_of_atoms = self.shx.atoms.number
        # This is a rough aproximation, but it works until you have a really high number of processors:
        barray = number_of_atoms * 8
        if barray <= 3000:
            barray = 3000
        return barray

    def remove_acta_card(self, acta_card):
        """
        Removes ACTA x from reslist and stores value in self._acta_card.
        """
        self._acta_card = acta_card._textline.strip('\r\n')[:]
        del self.shx._reslist[self.shx.index_of(acta_card)]
        self.shx.acta = None

    def restore_acta_card(self):
        """
        Place ACTA after UNIT
        """
        if not self._acta_card:
            return
        acta = ACTA(self.shx, self._acta_card.split())
        self.shx._reslist.insert(self.shx.unit.index + 1, ' ')
        self.shx.acta = self.shx.assign_card(acta, self.shx.unit.index + 1)

    def backup_shx_file(self):
        """
        makes a copy of the res file
        make backup in shxsaves before every fragment fit.
        name: self.resfile_name-date-time-seconds.res
        """
        import datetime
        now = datetime.datetime.now()
        timestamp = (str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '_' +
                     str(now.hour) + '-' + str(now.minute) + '-' + str(now.second))
        resfile = os.path.abspath(str(self.resfile_name + '.res'))
        bakup_dir = os.path.abspath(os.path.dirname(resfile)) + os.path.sep + os.path.relpath('shxsaves')
        try:
            copyfile(resfile, self.backup_file)
        except IOError:
            print('*** Unable to make backup file from {}. ***'.format(resfile))
            sys.exit()
        if not os.path.exists(bakup_dir):
            try:
                os.makedirs(bakup_dir)
            except(IOError, OSError):
                print('*** Unable to create backup directory {}. ***'.format(bakup_dir))
        try:
            copyfile(resfile,
                     bakup_dir + os.path.sep + os.path.split(self.resfile_name)[1] + '_' + timestamp + '.res')
        except IOError:
            print('\n*** Unable to make backup file from {} in shxsaves. ***'.format(resfile))

    def restore_shx_file(self):
        """
        restores filename from backup
        """
        resfile = os.path.abspath(str(self.resfile_name + '.res'))
        try:
            print('*** Restoring previous res file. ***')
            copyfile(self.backup_file, resfile)
        except IOError:
            print('Unable to restore res file from {}.'.format(self.backup_file))
        try:
            remove_file(self.backup_file)
        except IOError:
            print('Unable to delete backup file {}.'.format(self.backup_file))

    def pretty_shx_output(self, output):
        """
        selectively prints the output from shelx
        """
        for out in output:
            if out.startswith(' +  Copyright(C)'):
                print(' SHELXL {}'.format(' '.join(out.split()[6:8])))
            # wR2
            # These values are always bad after a simple LS fit without any atom movement:
            # if out.startswith(' wR2') and not wr2:
            #    wr2 = True
            #    line = out[:].split()
            #    print(' {}  {} {:>6}'.format(line[0], line[1], line[2][:6]))
            # R1
            # if out.startswith(' R1') and not r1:
            #    r1 = True
            #    line = out[:].split()
            #    print(' {}   {} {:>6}'.format(line[0], line[1], line[2][:6]))
            # GooF
            # if re.match(r'.*GooF.*', out) and not gof:
            #    gof = True
            #    line = out.split()
            #    print(' {} {} {:>5}0'.format(line[0], line[1], line[4][:5]))
            if re.match(r'.*CANNOT RESOLVE (SAME|RIGU|SIMU|DELU)', out):
                print('\nWarning: Are you sure that all atoms are in the correct order?\n')
            if re.match(r'.*CANNOT\s+OPEN\s+FILE.*hkl.*', out):
                print('*** No hkl file found! ***')
                print('*** You need a proper hkl file to run SHELXL! ***')
                sys.exit()
            if re.match(r'.*\*\* Extinction \(EXTI\) or solvent.*', out):
                continue
            if re.match(r'.*\*\* MERG code changed to 0', out):
                # disable this output
                continue
            if re.match(r'.*\*\* Bond\(s\) to .* ignored', out):
                # disable this output
                continue
            if re.match(r'.*\*\*.*', out):
                print('\n SHELXL says:')
                print(' {}'.format(out.strip('\n\r')))

    def run_shelxl(self):
        """
        This method runs shelxl 2013 on the res file self.resfile_name
        """
        resfile = self.resfile_name + '.res'
        hklfile = self.resfile_name + '.hkl'
        if not os.path.exists(hklfile):
            print('You need a proper hkl file to run SHELXL.')
            sys.exit()
        command_line = ['{}'.format(self._shelx_command), "-b{}".format(self.get_b_array()),
                        '{}'.format(self.resfile_name)]
        self.backup_shx_file()
        print(sep_line)
        print(' Running SHELXL with "{}" and "{}"'.format(' '.join(command_line), self.shx.cycles))
        p = subprocess.Popen(command_line,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)
        child_stdin.close()
        # Watch the output for successful termination
        out = child_stdout_and_stderr.readline().decode('ascii')
        output = []
        while out:
            if "cycle" in out:
                print(out.rstrip("\n\r"))
            output.append(out)
            out = child_stdout_and_stderr.readline().decode('ascii')
        child_stdout_and_stderr.close()
        # output only the most importand things from shelxl:
        self.pretty_shx_output(output)
        status = True
        if os.stat(resfile).st_size < 10:
            # status is False if shelx was unsecessful
            status = False
        if not status:  # fail
            print(sep_line)
            print('\nError: SHELXL terminated unexpectedly.')
            print('Check for errors in your SHELX input file!\n')
            self.restore_shx_file()
            sys.exit()

    def check_refinement_results(self, list_file):
        """
        Does some checks if the refinement makes sense e.g. if the data to parameter
        ratio is in an acceptable range.
        :param list_file: SHELXL listing file
        :type list_file: list
        """
        is_resfile_there = os.path.exists(self.resfile_name + '.res')
        if is_resfile_there and is_resfile_there == 'zero':
            print('Something failed in SHELXL. Please check your .ins and .lst file!')
            self.restore_shx_file()
            try:
                remove_file(self.backup_file)
            except IOError:
                print('Unable to delete backup file {}.'.format(self.backup_file))
            sys.exit()
        if not is_resfile_there:
            print('Something failed in SHELXL. Please check your .ins and .lst file!')
            self.restore_shx_file()
            try:
                remove_file(self.backup_file)
            except IOError:
                print('Unable to delete backup file {}.'.format(self.backup_file))
            sys.exit()
        regex_final = r' Final Structure Factor Calculation.*\n'
        final_results = find_line(list_file, regex_final)
        # find data and parameters:
        try:
            dataobj = re.search(r'[0-9]+\s+data', list_file[final_results + 4])
            data = float(dataobj.group(0).split()[0])
            parameterobj = re.search(r'[0-9]+\s+parameters', list_file[final_results + 4])
            parameters = float(parameterobj.group(0).split()[0])
            restrobj = re.search(r'[0-9]+\s+restraints', list_file[find_line(list_file, r" GooF = S =.*")])
            restraints = float(restrobj.group(0).split()[0])
        except AttributeError:
            return False
        try:
            data_to_parameter_ratio = data / parameters
            restr_ratio = ((data + restraints) / parameters)
        except ZeroDivisionError:
            return False
        lattline = find_line(list_file, r'^ LATT.*')
        centro = None
        if lattline:
            try:
                latt = int(list_file[lattline].split()[1])
            except ValueError:
                latt = 1
            if latt > 0:
                centro = True
            else:
                centro = False
        if centro and data_to_parameter_ratio < 10:
            print('*** Warning! The data/parameter ratio is getting low (ratio = {:.1f})! ***'
                  '\n*** but consider (data+restraints)/parameter = {:.1f} ***'
                  .format(data_to_parameter_ratio, restr_ratio))
        if not centro and data_to_parameter_ratio < 7.5:
            print('*** Warning! The data/parameter ratio is getting low (ratio = {:.1f})! ***'
                  '\n*** but consider (data+restraints)/parameter = {:.1f} ***'
                  .format(data_to_parameter_ratio, restr_ratio))
        try:
            remove_file(self.backup_file)
        except IOError:
            print('Unable to delete backup file {}.'.format(self.backup_file))


if __name__ == '__main__':
    pass
