import numpy as np
import matplotlib.pyplot as plt
import copy

class Keyword(object):
    def __init__(self, write_file_name):
        self.write_file_name = write_file_name
        self.lines = []
        self.lines.append(['*KEYWORD LONG=Y\n'])

    ########################################################
    # Blocks
    ########################################################

    def node(self, nodes, tc=0, rc=0):
        if nodes.__class__ == {}.__class__:
            nodes_iter=iter(nodes.values())
        elif nodes.__class__ == [].__class__:
            nodes_iter=nodes
        else:
            raise Exception('Node list not list or dict')
        line_block = ['*Node\n']
        line_block.append(self.format_comment_line(['nid', 'x', 'y', 'z', 'tc', 'rc']))
        for node in nodes_iter:
            line_block.append(self.format_key_line([node.id_, *list(node.coord), tc, rc]))
        self.submit_block(line_block)

    def element_shell(self, elements):
        if elements.__class__ == {}.__class__:
            elements_iter = iter(elements.values())
        elif elements.__class__ == [].__class__:
            elements_iter = elements
        else:
            raise Exception('Element list not list or dict')
        line_block = ['*ELEMENT_SHELL\n']
        line_block.append(self.format_comment_line(['nid', 'pid', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6']))
        for elem in elements_iter:
            line_block.append(self.format_key_line([elem.id_, elem.parent, *elem.node_ids]))
        self.submit_block(line_block)

    def element_solid(self, elements):
        if elements.__class__ == {}.__class__:
            elements_iter = iter(elements.values())
        elif elements.__class__ == [].__class__:
            elements_iter = elements
        else:
            raise Exception('Element list not list or dict')
        line_block = ['*ELEMENT_SOLID\n']
        line_block.append(self.format_comment_line(['nid', 'pid', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6']))
        for elem in elements_iter:
            line_block.append(self.format_key_line([elem.id_, elem.parent, *elem.node_ids]))
        self.submit_block(line_block)

    def element_shell_offset(self, elements, offset):
        if elements.__class__ == {}.__class__:
            elements_iter = iter(elements.values())
        elif elements.__class__ == [].__class__:
            elements_iter = elements
        else:
            raise Exception('Element list not list or dict')
        line_block = ['*ELEMENT_SHELL_OFFSET\n']
        line_block.append(self.format_comment_line(['nid', 'pid', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6']))
        for elem in elements_iter:
            line_block.append(self.format_key_line([elem.id_, elem.parent, *elem.node_ids]))
            line_block.append(self.format_key_line([offset]))
        self.submit_block(line_block)

    def element_beam_orientation(self, elements):
        if elements.__class__ == {}.__class__:
            elements_iter = iter(elements.values())
        elif elements.__class__ == [].__class__:
            elements_iter = elements
        else:
            raise Exception('Element list not list or dict')
        line_block = ['*ELEMENT_BEAM_ORIENTATION\n']
        line_block.append(self.format_comment_line(['nid', 'pid', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6']))
        for elem in elements_iter:
            line_block.append(self.format_key_line([elem.id_, elem.parent, *elem.node_ids]))
            line_block.append(self.format_key_line([*list(elem.orientation)]))
        self.submit_block(line_block)

    def element_beam_thickness_orientation(self, elements):
        if elements.__class__ == {}.__class__:
            elements_iter = iter(elements.values())
        elif elements.__class__ == [].__class__:
            elements_iter = elements
        else:
            raise Exception('Element list not list or dict')
        line_block = ['*ELEMENT_BEAM_THICKNESS_ORIENTATION\n']
        line_block.append(self.format_comment_line(['nid', 'pid', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6']))
        for elem in elements_iter:
            ts1 = np.sqrt(elem.csa * 4 / np.pi)
            tt1 = 0
            line_block.append(self.format_key_line([elem.id_, elem.parent, *elem.node_ids]))
            line_block.append(self.format_key_line([ts1, ts1, tt1, tt1,]))
            line_block.append(self.format_key_line([*list(elem.orientation)]))
        self.submit_block(line_block)

    def element_beam_section07_orientation(self, elements):
        if elements.__class__ == {}.__class__:
            elements_iter = iter(elements.values())
        elif elements.__class__ == [].__class__:
            elements_iter = elements
        else:
            raise Exception('Element list not list or dict')
        line_block = ['*ELEMENT_BEAM_SECTION_ORIENTATION\n']
        line_block.append(self.format_comment_line(['nid', 'pid', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6']))
        stype = 'SECTION_07'
        for elem in elements_iter:
            D1 = np.sqrt(elem.csa * 4 / np.sqrt(3))
            D2 = D1 / 100
            D3 = D1 * np.sqrt(3) / 2
            D4 = 0
            line_block.append(self.format_key_line([elem.id_, elem.parent, *elem.node_ids]))
            line_block.append(self.format_comment_line(['stype', 'D1', 'D2', 'D3', 'D4', 'D5']))
            line_block.append(self.format_key_line([stype, D1, D2, D3, D4]))
            line_block.append(self.format_comment_line([' vx', 'vy', 'vz']))
            line_block.append(self.format_key_line([*list(elem.orientation)]))
        self.submit_block(line_block)

    def end_key(self):
        self.submit_block(['*END\n'])

    ########################################################
    # Read and write
    ########################################################
    def write_key(self):
        flatt_lines = []
        for sublist in self.lines:
            for line in sublist:
                flatt_lines.append(line)
        with open(self.write_file_name, 'w+') as write_file:
            write_file.writelines(flatt_lines)

    ########################################################
    # Bookkeeping
    ########################################################
    def comment_block(self, commentStr, insertInd=False):
        line_block = ['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n']
        line_block.append('$')
        line_block.append('$$$$    {}\n'.format(commentStr))
        line_block.append('$')
        line_block.append('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        self.lines.append(line_block)

    def format_key_line(self, list_of_items):
        # comma separated, max 8 integers (I8)
        string_array = ''
        for item in list_of_items:
            if isinstance(item, float):
                if abs(item) < 1e6 and abs(item) > 1e-4:
                    string_array = string_array + '{:>20s}'.format(str(item)[:19])
                elif item == 0.0:
                    string_array = string_array + '{:>20.1f}'.format(item)
                else:
                    string_array = string_array + '{:>20.14g}'.format(item)
            if isinstance(item, int):
                string_array = string_array + '{:>20.0f}'.format(item)
            if isinstance(item, str):
                string_array = string_array + '{:<20s}'.format(item)
        return string_array + '\n'

    def format_comment_line(self, list_of_items):
        # comma separated, max 8 integers (I8)
        string_array = '$#'+list_of_items[0].rjust(18, ' ')
        for item in list_of_items[1:]:
            string_array = string_array + item.rjust(20, ' ')
        return string_array + '\n'

    def format_key_line_short(self, list_of_items):
        # comma separated, max 8 integers (I8)
        string_array = ''
        for item in list_of_items:
            if isinstance(item, float):
                if abs(item) < 1e6 and abs(item) > 1e-4:
                    string_array = string_array + '{:>10s}'.format(str(item)[:9])
                elif item == 0.0:
                    string_array = string_array + '{:>10.1f}'.format(item)
                else:
                    string_array = string_array + '{:>10.14g}'.format(item)
            if isinstance(item, int):
                string_array = string_array + '{:>10.0f}'.format(item)
            if isinstance(item, str):
                string_array = string_array + '{:<10s}'.format(item)
        return string_array + '\n'

    def submit_block(self, line_block):
        self.lines.append(line_block)

    def set_node_list(self, nsid, node_list):
        line_block = ['*SET_NODE_LIST_TITLE\n']
        line_block.append('Node list {}\n'.format(nsid))
        line_block.append(self.format_comment_line(['sid', 'da1', 'da2', 'da3', 'da4', 'solver']))
        line_block.append(self.format_key_line([nsid, 0.0, 0.0, 0.0, 0.0, 'MECH']))
        line_block.append(self.format_comment_line(['nid1', 'nid2', 'nid3', 'nid4', 'nid5', 'nid6', 'nid7', 'nid8']))
        node_list=copy.copy(node_list)
        k, m = divmod(len(node_list), 8)
        if m != 0:
            node_list.extend([0] * (8 - m))
            k = k + 1
        for i in range(k):
            line_block.append(self.format_key_line(node_list[8 * i:8 * (i + 1)]))
        self.submit_block(line_block)

    def set_segment(self, nsid, node_list):
        line_block = ['*SET_SEGMENT\n']
        line_block.append(self.format_comment_line(['sid', 'da1', 'da2', 'da3', 'da4', 'solver']))
        line_block.append(self.format_key_line([nsid, 0.0, 0.0, 0.0, 0.0, 'MECH']))
        line_block.append(self.format_comment_line(['nid1', 'nid2', 'nid3', 'nid4', 'nid5', 'nid6', 'nid7', 'nid8']))
        node_list = copy.copy(node_list)
        k, m = divmod(len(node_list), 8)
        if m != 0:
            node_list.extend([0] * (8 - m))
            k = k + 1
        for i in range(k):
            line_block.append(self.format_key_line(node_list[8 * i:8 * (i + 1)]))
        self.submit_block(line_block)

    def set_2D_segment(self, sid, pid_list):
        line_block = ['*SET_2D_SEGMENT_TITLE\n']
        line_block.append('Volume set {}\n'.format(sid))
        line_block.append(self.format_comment_line(['sid', 'da1', 'da2', 'da3', 'da4']))
        line_block.append(self.format_key_line([sid, 0.0, 0.0, 0.0, 0.0]))
        line_block.append(self.format_comment_line(['pid']))
        for pid in pid_list:
            line_block.append(self.format_key_line([pid]))
        self.submit_block(line_block)

    def set_part_list(self, sid, pid_list):
        line_block = ['*SET_PART_LIST\n']
        line_block.append(self.format_comment_line(['sid', 'da1', 'da2', 'da3', 'da4', 'solver']))
        line_block.append(self.format_key_line([sid, 0.0, 0.0, 0.0, 0.0, 'MECH']))
        line_block.append(self.format_comment_line(['pid1', 'pid2', 'pid3', 'pid4', 'pid5', 'pid6', 'pid7', 'pid8']))
        pid_list=copy.copy(pid_list)
        k, m = divmod(len(pid_list), 8)
        if m != 0:
            pid_list.extend([0] * (8 - m))
            k = k + 1
        for i in range(k):
            line_block.append(self.format_key_line(pid_list[8 * i:8 * (i + 1)]))
        self.submit_block(line_block)

    ###########################################
    # Control
    ################################################
    def control_termination(self, endtim, dtmin=0.0):
        self.endtim = endtim
        endcyc = 0
        endeng = 0.0
        endmas = 1.0E08  # E08
        nosol = 0
        line_block = ['*CONTROL_TERMINATION\n']
        line_block.append(self.format_comment_line(['endtim', 'endcyc', 'dtmin', 'endeng', 'endmas', 'nosol']))
        line_block.append(self.format_key_line([endtim, endcyc, dtmin, endeng, endmas, nosol]))
        self.submit_block(line_block)

    def control_structured(self):
        line_block = ['*CONTROL_STRUCTURED\n']
        self.submit_block(line_block)

    def control_accuracy(self, osu=1, inn=2, iacc=1):
        line_block = ['*CONTROL_ACCURACY\n']
        line_block.append(self.format_comment_line(['osu', 'inn', 'pidosu', 'iacc']))
        line_block.append(self.format_key_line([osu, inn, 0, iacc]))
        self.submit_block(line_block)

    def control_implicit_dynamics(self, imass=1, gamma=0.6, beta=0.38):
        line_block = ['*CONTROL_IMPLICIT_DYNAMICS\n']
        line_block.append(self.format_comment_line(['imass', 'gamma', 'beta', 'tdybir', 'tdydth', 'tdybur', 'irate', 'alpha']))
        line_block.append(self.format_key_line([imass, gamma, beta, 0.0, 9e27, 9e27, 0, 0.0]))
        self.submit_block(line_block)

    def control_implicit_solution(self, nsolvr=12, ilimit=11, maxref=15, dctol=0.001, ectol=0.01):
        line_block = ['*CONTROL_IMPLICIT_SOLUTION\n']
        line_block.append(self.format_comment_line(['nsolvr', 'ilimit', 'maxref', 'dctol', 'ectol', 'rctol', 'lstol', 'abstol']))
        line_block.append(self.format_key_line([nsolvr, ilimit, maxref, dctol, ectol, 1e10, 0.9, 1e-10]))
        line_block.append(self.format_comment_line(['dnorm', 'diverg', 'istif', 'nlprint', 'nlnorm', 'd3itctl', 'cpchk']))
        line_block.append(self.format_key_line([2, 1, 1, 0, 2, 0, 0]))
        line_block.append(self.format_comment_line(['arcctl', 'arcdir', 'arclen', 'arcmth', 'arcdmp', 'arcpsi', 'arcalf', 'arctim']))
        line_block.append(self.format_key_line([0, 0, 0.0, 1, 2, 0, 0, 0]))
        line_block.append(self.format_comment_line(['lsmtd', 'lsdir', 'irad', 'srad', 'awgt', 'sred']))
        line_block.append(self.format_key_line([4, 2, 0.0, 0.0, 0.0, 0.0]))
        self.submit_block(line_block)

    def control_implicit_solver(self):
        line_block = ['*CONTROL_IMPLICIT_SOLVER\n']
        line_block.append(self.format_comment_line(['lsolvr', 'lprint', 'negev', 'order', 'drcm', 'drcprm', 'autospc', 'autotol']))
        line_block.append(self.format_key_line([2, 0, 2, 0, 4, 0.0, 1, 0.0]))
        line_block.append(self.format_comment_line(['lcpack', 'mtxdmp', 'iparm1', 'rparm1', 'rparm2']))
        line_block.append(self.format_key_line([2, 0, 500, 1e-9, 1e-3]))
        line_block.append(self.format_comment_line(['emxdmp', 'rdcmem']))
        line_block.append(self.format_key_line([0, 0.85]))
        self.submit_block(line_block)

    def control_shell(self, istupd=4, psstupd=0, irnxx=-2, miter=2, esort=1, nfail1=1, nfail4=1):
        '''esort 1 sorts degenerate quadrialateral elements as triangular C0 elements'''
        wrpang = 20.0

        intgrd = 0  # Gauss integration
        cntco = 1

        line_block = ['*CONTROL_SHELL\n']
        line_block.append(self.format_comment_line(['wrpang', 'esort', 'irnxx', 'istupd', 'theory', 'bwc', 'miter', 'proj']))
        line_block.append(self.format_key_line([wrpang, esort, irnxx, istupd, 2, 1, miter, 1]))
        line_block.append(self.format_comment_line(['rotascl', 'intgrd', 'lamsht', 'cstyp6', 'thshel']))
        line_block.append(self.format_key_line([1.0, intgrd, 0, 1, 0]))
        line_block.append(self.format_comment_line(['psstupd', 'sidt4tu', 'cntco', 'itsflg', 'irquad', 'w-mode', 'stretch', 'icrq']))
        line_block.append(self.format_key_line([psstupd, 0, cntco, 0, 2, 0.0000000000e+000, 0.0000000000e+000, 0]))
        line_block.append(self.format_comment_line(['nfail1', 'nfail4', 'psnfail', 'keepcs', 'delfr', 'drcpsid', 'drcprm', 'intperr']))
        line_block.append(self.format_key_line([nfail1, nfail4, 0, 0, 0, 0, 1.0, 0]))
        self.submit_block(line_block)

    def control_timestep(self, dt2ms, tssfac=0.9):
        '''TSSFAC*dt2ms is the minimun allowed timestep'''

        line_block = ['*CONTROL_TIMESTEP\n']
        line_block.append(self.format_comment_line(['dtinit', 'tssfac', 'isdo', 'tslimt', 'dt2ms', 'lctm', 'erode', 'ms1st']))
        line_block.append(self.format_key_line([0.0, tssfac, 0, 0.0, dt2ms, 0, 0, 0]))
        line_block.append(self.format_comment_line(['dt2msf', 'dt2mslc', 'imscl', 'unused', 'unused', 'rmscl']))
        line_block.append(self.format_key_line([0.0, 0, 0, '', '', 0.0]))
        self.submit_block(line_block)

    def control_contact(self, shlthk=2):
        slsfac = 0.1
        rwpnal = 0.0
        islchk = 1
        penopt = 1
        thkchg = 1
        orien = 1
        enmass = 0
        ssthk = 1
        usrstr, usrfrc, nsbcs, interm, xpene, ecdt, tiedprj = 0, 0, 0, 0, 4.0, 0, 0
        line_block = ['*CONTROL_CONTACT\n']
        line_block.append(self.format_comment_line(['slsfac', 'rwpnal', 'islchk', 'shlthk', 'penopt', 'thkchg', 'orien', 'enmass']))
        line_block.append(self.format_key_line([slsfac, rwpnal, islchk, shlthk, penopt, thkchg, orien, enmass]))
        line_block.append(self.format_comment_line(['usrstr', 'usrfrc', 'nsbcs', 'interm', 'xpene', 'ssthk', 'ecdt', 'tiedprj']))
        line_block.append(self.format_key_line([usrstr, usrfrc, nsbcs, interm, xpene, ssthk, ecdt, tiedprj]))
        line_block.append(self.format_comment_line(['sfric', 'dfric', 'edc', 'vfc', 'th', 'th_sf', 'pen_sf']))
        line_block.append(self.format_key_line([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        line_block.append(self.format_comment_line(['ignore', 'frceng', 'skiprwg', 'outseg', 'spotstp', 'spotdel', 'spothin']))
        line_block.append(self.format_key_line([1, 0, 0, 0, 0, 0, 0.0]))
        line_block.append(self.format_comment_line(['isym', 'nserod', 'rwgaps', 'rwgdth', 'rwksf', 'icov', 'swradf', 'ithoff']))
        line_block.append(self.format_key_line([0, 1, 1, 0.0, 1.0, 0, 0.0, 1]))
        line_block.append(self.format_comment_line(['shledg', 'pstiff', 'ithcnt', 'tdcnof', 'ftall', 'unused', 'shltrw']))
        line_block.append(self.format_key_line([0, 0, 0, 1, 0, 0, 0.0]))
        self.submit_block(line_block)

    def control_hourglass(self, ihg=4, qh=0.1):

        line_block = ['*CONTROL_HOURGLASS\n']
        line_block.append(self.format_comment_line(['ihq', 'qh']))
        line_block.append(self.format_key_line([ihg, qh]))
        self.submit_block(line_block)

    def control_implicit_general(self, imflag, dt0=None):
        if dt0 == None:
            dt0 = self.endtim / 10.
        line_block = ['*CONTROL_IMPLICIT_GENERAL\n']
        line_block.append(self.format_comment_line(['imflag', 'dt0', 'imform', 'nsbs', 'igs', 'cnstn', 'form', 'zero_v']))
        line_block.append(self.format_key_line([imflag, dt0, 2, 1, 2, 0, 0, 0]))
        self.submit_block(line_block)

    def control_implicit_auto(self, iauto=1, dtmin=0.0, dtmax=0.0, iteopt=11, dtexp=1e-9):
        itewin = 5
        line_block = ['*CONTROL_IMPLICIT_AUTO\n']
        line_block.append(self.format_comment_line(['iauto', 'iteopt', 'itewin', 'dtmin', 'dtmax', 'dtexp', 'kfail', 'kcycle']))
        line_block.append(self.format_key_line([iauto, iteopt, itewin, dtmin, dtmax, dtexp, 0, 0]))
        self.submit_block(line_block)

    def control_parallel(self, ncpu, const=1):
        line_block = ['*CONTROL_PARALLEL\n']
        line_block.append(self.format_comment_line(['ncpu', 'numrhs', 'const', 'para']))
        line_block.append(self.format_key_line([ncpu, 0, const, 1]))
        self.submit_block(line_block)

    ###########################################
    # Database
    ################################################
    def database_glstat(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_GLSTAT\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_binary_d3_plot(self, dt=None):
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_BINARY_D3PLOT\n']
        line_block.append(self.format_comment_line(['dt', 'lcdt', 'beam', 'psetid']))
        line_block.append(self.format_key_line([dt, 0, 0, 0]))
        line_block.append(self.format_comment_line(['ioopt']))
        line_block.append(self.format_key_line([ioopt]))
        self.submit_block(line_block)

    def database_rcforc(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_RCFORC\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_abstat(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_ABSTAT\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_rwforc(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_RWFORC\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_rbdout(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_RBDOUT\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_spcforc(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_SPCFORC\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_bndout(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_BNDOUT\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_nodout(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        option1 = 0
        option2 = 0
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_NODOUT\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt', 'option1', 'option2']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt, option1, option2]))
        self.submit_block(line_block)

    def database_nodfor(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_NODFOR\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_nodfor_group(self, nsid):
        line_block = ['*DATABASE_NODAL_FORCE_GROUP\n']
        line_block.append(self.format_comment_line(['nsid', 'cid']))
        line_block.append(self.format_key_line([nsid, 0]))
        self.submit_block(line_block)

    def database_ncforc(self, dt=None):
        binary = 1
        lcur = 0
        ioopt = 1
        if dt == None:
            dt = self.endtim / 25.
        line_block = ['*DATABASE_NCFORC\n']
        line_block.append(self.format_comment_line(['dt', 'binary', 'lcur', 'ioopt']))
        line_block.append(self.format_key_line([dt, binary, lcur, ioopt]))
        self.submit_block(line_block)

    def database_hist_node_set(self, nsid):
        line_block = ['*DATABASE_HISTORY_NODE_SET\n']
        line_block.append(self.format_comment_line(['id1', 'id2', 'id3', 'id4', 'id5', 'id6', 'id7', 'id8']))
        line_block.append(self.format_key_line(nsid))
        self.submit_block(line_block)

    def database_hist_node(self, nids):
        line_block = ['*DATABASE_HISTORY_NODE\n']
        line_block.append(self.format_comment_line(['id1', 'id2', 'id3', 'id4', 'id5', 'id6', 'id7', 'id8']))
        nids = copy.copy(nids)
        k, m = divmod(len(nids), 8)
        if m != 0:
            nids.extend([0] * (8 - m))
            k = k + 1
        for i in range(k):
            line_block.append(self.format_key_line(nids[8 * i:8 * (i + 1)]))
        self.submit_block(line_block)

    ###########################################
    # Parts, sections
    ################################################

    def part(self, pid, secid, mid, eosid=0):
        hgid = 0
        gray = 0
        adpopt = 0
        tmid = 0
        line_block = ['*PART\n']
        line_block.append('Auto from elem')
        line_block.append(self.format_comment_line(['pid', 'secid', 'mid', 'eosid', 'hgid', 'grav', 'adpopt', 'tmid']))
        line_block.append(self.format_key_line([pid, secid, mid, eosid, hgid, gray, adpopt, tmid]))
        self.submit_block(line_block)

    def part_contact(self, pid, secid, mid, eosid=0, fs=0.0, fd=0.0, dc=0.0):
        hgid = 0
        gray = 0
        adpopt = 0
        tmid = 0
        line_block = ['*PART_CONTACT\n']
        line_block.append('Contact')
        line_block.append(self.format_comment_line(['pid', 'secid', 'mid', 'eosid', 'hgid', 'grav', 'adpopt', 'tmid']))
        line_block.append(self.format_key_line([pid, secid, mid, eosid, hgid, gray, adpopt, tmid]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'opt', 'sft', 'ssf', 'CPARM8']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        self.submit_block(line_block)

    def section_shell(self, secid, t1, elform=1, nip=5, nloc=0.0):
        shrf = 1.0
        propt = 1.0
        qririd = 0
        icomp = 0
        setyp = 1
        marea = 0.0
        idof = 0.0
        edgset = 0
        line_block = ['*SECTION_SHELL_TITLE\n']
        line_block.append('Section {}\n'.format(secid))
        line_block.append(self.format_comment_line(['secid', 'elform', 'shrf', 'nip', 'propt', 'qr/irid', 'icomp', 'setyp']))
        line_block.append(self.format_key_line([secid, elform, shrf, nip, propt, qririd, icomp, setyp]))

        line_block.append(self.format_comment_line(['t1', 't2', 't3', 't4', 'nloc', 'marea', 'idof', 'edgset']))
        line_block.append(self.format_key_line([t1, t1, t1, t1, nloc, marea, idof, edgset]))

        self.submit_block(line_block)

    def section_solid(self, secid, elform=1):
        line_block = ['*SECTION_SOLID_TITLE\n']
        line_block.append('Section {}\n'.format(secid))
        line_block.append(self.format_comment_line(['secid', 'elform', 'aet']))
        line_block.append(self.format_key_line([secid, elform, 0]))
        self.submit_block(line_block)

    def section_beam(self, secid, csa=1.0, elform=1):
        shrf = 1.0
        qririd = 5.0
        cst = 1.0
        scoor = 0.0
        nsm = 0.0
        ts1 = np.sqrt(csa * 4 / np.pi)
        tt1 = 0
        nsloc = 0.0
        ntloc = 0.0
        line_block = ['*SECTION_BEAM_TITLE\n']
        line_block.append('Section {}\n'.format(secid))
        line_block.append(self.format_comment_line(['secid', 'elform', 'shrf', 'QR/IRID', 'cst', 'scoor', 'nsm']))
        line_block.append(self.format_key_line([secid, elform, shrf, qririd, cst, scoor, nsm]))

        line_block.append(self.format_comment_line(['ts1', 'ts2', 'tt1', 'tt2', 'nsloc', 'ntloc']))
        line_block.append(self.format_key_line([ts1, ts1, tt1, tt1, nsloc, ntloc]))

        self.submit_block(line_block)
    ###########################################
    # Material
    ################################################
    def mat24(self, mid, ro=9.20000E-10, e=1500.0, pr=0.3, sigy=20.,
                               etan=1., fail=0.0, tdel=0.0, lcss=0, c=0, p=0):

        line_block = ['*MAT_PIECEWISE_LINEAR_PLASTICITY_TITLE\n']
        line_block.append('linear - plasticity')
        line_block.append(self.format_comment_line(['mid', 'ro', 'e', 'pr', 'sigy', 'etan', 'fail', 'tdel']))
        line_block.append(self.format_key_line([mid, ro, e, pr, sigy, etan, fail, tdel]))

        line_block.append(self.format_comment_line(['c', 'p', 'lcss', 'lcsr', 'vp']))
        line_block.append(self.format_key_line([c, p, lcss, 0, 0.0]))

        line_block.append(self.format_comment_line(['eps1', 'eps2', 'eps3', 'eps4', 'eps5', 'eps6', 'eps7', 'eps8']))
        line_block.append(self.format_key_line([0.0] * 8))

        line_block.append(self.format_comment_line(['es1', 'es2', 'es3', 'es4', 'es5', 'es6', 'es7', 'es8']))
        line_block.append(self.format_key_line([0.0] * 8))

        self.submit_block(line_block)

    def mat28(self, mid, ro=9.20000E-10, e=1500.0, pr=0.3, sigy=20., etan=1.):
        line_block = ['*MAT_RESULTANT_PLASTICITY_TITLE\n']
        line_block.append('linear - plasticity')
        line_block.append(self.format_comment_line(['mid', 'ro', 'e', 'pr', 'sigy', 'etan', 'fail', 'tdel']))
        line_block.append(self.format_key_line([mid, ro, e, pr, sigy, etan]))
        self.submit_block(line_block)

    def mat24_rate(self, mid, ro=9.20000E-10, e=1500.0, pr=0.3, sigy=20.,
                        etan=1., fail=0.0, tdel=0.0, str_mod=[0.075, 0.2, 1e-3, 3e1], plot_curve=False, soften=False):
        tbid = 2000
        sigma_y = sigy / (1 + str_mod[0] * np.log10(str_mod[2]))

        strain_rate = np.logspace(np.log10(str_mod[2]), np.log10(str_mod[2]) + 7, 20)

        a = 1 + str_mod[0] * np.log10(strain_rate)
        # 1+str_mod[0] * np.log10(strain_rate) = b + str_mod[1] * np.log10(strain_rate)
        # b =  1 + (str_mod[0] * np.log10(str_mod[3])) - (str_mod[1] * np.log10(str_mod[3]))
        b = 1 + (str_mod[0] - str_mod[1]) * np.log10(str_mod[3])
        c = b + str_mod[1] * np.log10(strain_rate)
        scaled_yield_stress = sigma_y * np.maximum(a, c)

        if plot_curve == True:
            plt.figure()
            plt.xscale('log')
            plt.ylim(0, 120)
            plt.scatter(strain_rate, a * sigma_y)
            plt.scatter(strain_rate, c * sigma_y)
            plt.grid()

        curve_id_list = np.arange(2101, 2101 + len(strain_rate), 1)
        self.define_table(tbid=tbid, value_list=strain_rate, id_list=curve_id_list)

        for yield_stress, curve_id in zip(scaled_yield_stress, curve_id_list):
            if fail == 0.0:
                failCurve = 20
            else:
                failCurve = fail
            abscissa = np.linspace(0, failCurve + 1, 10)
            ordinate = yield_stress + etan * abscissa
            if soften  == True:
                mat_dict=self.soften_yield_curve(youngs=e, sigy=yield_stress, etan=etan)
                self.define_curve(lcid=int(curve_id), abscissas=mat_dict['strain'], ordinates=mat_dict['stress'])
            else:
                self.define_curve(lcid=int(curve_id), abscissas=abscissa, ordinates=ordinate)

        self.mat24(mid, ro=ro, e=e, pr=pr, sigy=sigy,
                       etan=etan, fail=fail, tdel=tdel, lcss=tbid, c=0, p=0)

    def soften_yield_curve(self, youngs=1500.0, sigy=20., etan=1., plot_curve=False):
        yieldStrain = sigy / youngs
        E = youngs
        sigma = np.linspace(0, sigy * 1.5, 200)
        eps1 = 0.5
        eps2 = 5.0
        n = ((sigy + eps2 * etan) / (sigy + eps1 * etan)) / (eps2 / eps1)
        K = (sigy + eps1 * etan) / (eps1 ** (n))

        eps = sigma / E + (sigma / K) ** (1 / n)
        epsPl = (sigma / K) ** (1 / n)

        if plot_curve == True:
            plt.figure()
            plt.plot([0, yieldStrain], [0, sigy], [yieldStrain, (sigy + etan * 3) / youngs + 3],
                     [sigy, sigy + etan * 3],
                     color='black')
            plt.plot(eps, sigma)
            plt.plot(epsPl, sigma)
            # plt.plot(matDict['strain'], matDict['stress'])

        mat_dict = {}
        mat_dict['e'] = E
        mat_dict['stress'] = np.interp(np.linspace(eps[(sigma > sigy / 2).argmax()], 10, 300), epsPl, sigma)
        mat_dict['strain'] = np.linspace(0, 10 - eps[(sigma > sigy / 2).argmax()], 300)
        return mat_dict

    def mat53(self, mid, lcid, ro=3.0E-11, e=20, a=0, b=0, c=0, p0=0.1, phi=0.0333, gama0=0.0):
        line_block = ['*MAT_CLOSED_CELL_FOAM\n']
        line_block.append(self.format_comment_line(['mid', 'ro', 'e', 'a', 'b', 'c', 'p0', 'phi']))
        line_block.append(self.format_key_line([mid, ro, e, a, b, c, p0, phi]))
        line_block.append(self.format_comment_line(['gama0', 'lcid']))
        line_block.append(self.format_key_line([gama0, lcid]))
        self.submit_block(line_block)

    def mat57(self, mid, lcid, ro=3.0E-11, e=20):
        tc = 1e20
        hu = 1.0
        beta = 0.
        damp = 0.
        shape = 1.0
        fail = 1.0
        bvflag = 0.0
        ed = 0.0
        betal = 0.0
        kcon = 0.0
        ref = 0.0

        line_block = ['*MAT_LOW_DENSITY_FOAM\n']
        line_block.append(self.format_comment_line(['mid', 'ro', 'e', 'lcid', 'tc', 'hu', 'beta', 'damp']))
        line_block.append(self.format_key_line([mid, ro, e, lcid, tc, hu, beta, damp]))
        line_block.append(self.format_comment_line(['shape', 'fail', 'bvflag', 'ed', 'beta1', 'kcon', 'ref']))
        line_block.append(self.format_key_line([shape, fail, bvflag, ed, betal, kcon, ref]))
        self.submit_block(line_block)

    def mat181(self, mid, lcid, youngs=1500, ro=9.2E-10):
        pr = 0.3
        k = youngs / (3 * (1 - 2 * pr))
        mu = 0.00
        g = 0
        sigf = 0
        tension = -1
        rtype = 0.0
        avgopt = 0.0
        sgl = 1.0
        sw = 1.0
        st = 1.0
        line_block = ['*MAT_SIMPLIFIED_RUBBER/FOAM_WITH_FAILURE_TITLE\n']
        line_block.append('Rubber')
        line_block.append(self.format_comment_line(['mid', 'ro', 'k', 'mu', 'g', 'sigf']))
        line_block.append(self.format_key_line([mid, ro, k, mu, g, sigf, 0.0, 0.0]))
        line_block.append(self.format_comment_line(['sgl', 'sw', 'st', 'lc/tbid', 'tension', 'rtype', 'avgopt']))
        line_block.append(self.format_key_line([sgl, sw, st, lcid, tension, rtype, avgopt, pr]))
        line_block.append(self.format_comment_line(['k', 'gama1', 'gama2', 'eh']))
        line_block.append(self.format_key_line([0., 0., 0., 0.]))
        line_block.append(self.format_comment_line(['lcunld', 'hu', 'shape', 'stol', 'visco']))
        line_block.append(self.format_key_line([0, 1.0, 0., 0., 0.]))
        self.submit_block(line_block)

    def mat63(self, mid, lcid, ro=3.0E-11, e=20, pr=0.0, tsc=1e20):
        damp = 0.0

        line_block = ['*MAT_CRUSHABLE_FOAM\n']
        line_block.append(self.format_comment_line(['mid', 'ro', 'e', 'pr', 'lcid', 'tsc', 'damp']))
        line_block.append(self.format_key_line([mid, ro, e, pr, lcid, tsc, damp]))
        self.submit_block(line_block)

    def mat83(self, mid, tbid, ro=3.0E-11, e=20, tc=1e20):
        kcon = 0.0
        fail = 1.0
        damp = 0.0
        line_block = ['*MAT_FU_CHANG_FOAM_TITLE\n']
        line_block.append('Fu chang foam')
        line_block.append(self.format_comment_line(['mid', 'ro', 'e', 'kcon', 'tc', 'fail', 'damp', 'tbid']))
        line_block.append(self.format_key_line([mid, ro, e, kcon, tc, fail, damp, tbid]))
        line_block.append(self.format_comment_line(['bvflag', 'sflag', 'rflag', 'tflag', 'pvid', 'sraf', 'ref', 'hu']))
        line_block.append(self.format_key_line([0., 1.0, 0.0, 0.0, 0., 0.0, 0.0, 0.0]))
        line_block.append(self.format_comment_line(['d0', 'n0', 'n1', 'n2', 'n3', 'c0', 'c1', 'c2']))
        line_block.append(self.format_key_line([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        line_block.append(self.format_comment_line(['c3', 'c4', 'c5', 'aij', 'sij', 'minr', 'maxr', 'shape']))
        line_block.append(self.format_key_line([0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0]))
        line_block.append(self.format_comment_line(['expon', 'riuld']))
        line_block.append(self.format_key_line([1.0, 0.0]))
        self.submit_block(line_block)

    def mat154(self, mid, ro=3.0E-11, e=20, pr=0.2, alpha=0, gamma=0, epsd=0, alpha2=0, beta=0, sigp=0, derfi=0,
               cfail=0):
        line_block = ['*MAT_DESHPANDE_FLECK_FOAM_TITLE\n']
        line_block.append('Deshpande Fleck foam')
        line_block.append(self.format_comment_line([' mid', 'ro', 'e', 'pr', 'alpha', 'gamma']))
        line_block.append(self.format_key_line([mid, ro, e, pr, alpha, gamma]))
        line_block.append(self.format_comment_line(['epsd', 'alpha2', 'beta', 'sigp', 'derfi', 'cfail']))
        line_block.append(self.format_key_line([epsd, alpha2, beta, sigp, derfi, cfail]))
        self.submit_block(line_block)

    def mat_null(self, mid, e=1500, ro=9.4600E-10):
        pc = 0.0
        mu = 0.0
        terod = 0.0
        cerod = 0.0
        pr = 0.0
        line_block = ['*MAT_NULL_TITLE\n']
        line_block.append('Outer face material')
        line_block.append(self.format_comment_line(['mid', 'ro', 'pc', 'mu', 'terod', 'cerod', 'ym', 'pr']))
        line_block.append(self.format_key_line([mid, ro, pc, mu, terod, cerod, e, pr]))
        self.submit_block(line_block)

    def mat_elastic(self, mid, e=1500, ro=9.4600E-10, pr=0.3):
        line_block = ['*MAT_ELASTIC_TITLE\n']
        line_block.append('Ref element material')
        line_block.append(self.format_comment_line(['mid', 'ro', 'e', 'pr']))
        line_block.append(self.format_key_line([mid, ro, e, pr]))
        self.submit_block(line_block)

    def mat_rigid(self, mid, ro=9.4600E-10, e=1500):
        pr = 0.3
        line_block = ['*MAT_RIGID_TITLE\n']
        line_block.append('Outer face material')
        line_block.append(self.format_comment_line(['mid', 'ro', 'E', 'pr']))
        line_block.append(self.format_key_line([mid, ro, e, pr, 0, 0, 0]))
        line_block.append(self.format_comment_line(['cmo', 'con1', 'con2']))
        line_block.append(self.format_key_line([0, 0, 0]))
        line_block.append(self.format_comment_line(['div']))
        line_block.append(self.format_key_line([0, 0, 0, 0, 0, 0]))
        self.submit_block(line_block)

    ###########################################
    # Coordinate systems
    ################################################
    def define_coordinate_nodes(self, cid, nodes, flag):
        '''n1: origo, n2: x-axis, n3 : in xy plane'''
        dir = 'X'
        line_block = ['*DEFINE_COORDINATE_NODES_TITLE\n']
        line_block.append('Local Coord')
        line_block.append(self.format_comment_line(['cid', 'n1', 'n2', 'n3', 'flag', 'dir']))
        line_block.append(self.format_key_line([cid, nodes[0], nodes[1], nodes[2], flag, dir]))
        self.submit_block(line_block)

    def define_curve(self, lcid, abscissas, ordinates, sfo=1.0, sfa=1.0, offo=0.0, offa=0.0, dattyp=0, lcint=0):
        '''shape = 'vel' accelerates in endTime/100'''
        # offa offsets in time, to the right

        sidr = 0

        line_block = ['*DEFINE_CURVE_TITLE\n']
        line_block.append('define curve {}\n'.format(lcid))
        line_block.append(self.format_comment_line(['lcid', 'sidr', 'sfa', 'sfo', 'offa', 'offo', 'dattyp', 'lcint']))
        line_block.append(self.format_key_line([lcid, sidr, sfa, sfo, offa, offo, dattyp, lcint]))
        line_block.append(self.format_comment_line(['a1', 'o1']))
        for abscissa, ordinate in zip(abscissas, ordinates):
            line_block.append(self.format_key_line([float(abscissa), ordinate]))

        self.submit_block(line_block)

    def define_curve_smooth(self, lcid, dist, tstart, tend, trise_frac=0.15, v0=0.0):
        trise = (tend - tstart) * trise_frac
        sidr = 0
        line_block = ['*DEFINE_CURVE_SMOOTH\n']
        line_block.append(self.format_comment_line(['lcid', 'sidr', 'dist', 'tstart', 'tend', 'trise', 'v0']))
        line_block.append(self.format_key_line([lcid, sidr, dist, tstart, tend, trise, v0]))
        self.submit_block(line_block)

    def define_table_2D(self, tbid, value_list, id_list, sfa=1.0, offa=0.0):
        '''shape = 'vel' accelerates in endTime/100'''
        # offa offsets in time, to the right

        line_block = ['*DEFINE_TABLE_2D\n']
        line_block.append(self.format_comment_line(['tbid', 'sfa', 'offa']))
        line_block.append(self.format_key_line([tbid, sfa, offa]))
        line_block.append(self.format_comment_line(['value', 'curve id']))
        for abscissa, ordinate in zip(value_list, id_list):
            line_block.append(self.format_key_line([float(abscissa), int(ordinate)]))

        self.submit_block(line_block)

    def define_table(self, tbid, value_list, id_list, sfa=1.0, offa=0.0):
        '''shape = 'vel' accelerates in endTime/100'''
        # offa offsets in time, to the right

        line_block = ['*DEFINE_TABLE\n']
        # line_block.append('define table {}\n'.format(tbid))
        line_block.append(self.format_comment_line(['tbid', 'sfa', 'offa']))
        line_block.append(self.format_key_line([tbid, sfa, offa]))
        line_block.append(self.format_comment_line(['value', 'curve id']))
        for abscissa, ordinate in zip(value_list, id_list):
            line_block.append(self.format_key_line([abscissa, int(ordinate)]))

        self.submit_block(line_block)

    ###########################################
    # Curves
    ###############################################
    def loadCurveDispToVelMultiple(self, saturationAmpl, endTime):
        startTime = 0.0
        rampFraction = 0.1
        loadCurveDict = {}
        loadCurveDict['time'] = np.array([0.0])
        loadCurveDict['velAmpl'] = np.array([0.0])
        loadCurveDict['disp'] = np.array([0.0])
        steps = len(saturationAmpl)
        timeStep = endTime / steps
        for i, ampl in enumerate(saturationAmpl):
            velAmpl = ampl / (timeStep * (1 - rampFraction))
            loadCurveDict['time'] = np.append(loadCurveDict['time'],
                                              [loadCurveDict['time'][-1] + timeStep * rampFraction,
                                               loadCurveDict['time'][-1] + timeStep * (1 - rampFraction),
                                               loadCurveDict['time'][-1] + timeStep])
            loadCurveDict['velAmpl'] = np.append(loadCurveDict['velAmpl'], [velAmpl, velAmpl, 0.0])
        loadCurveDict['time'] = np.append(loadCurveDict['time'], [endTime + 1.])
        loadCurveDict['velAmpl'] = np.append(loadCurveDict['velAmpl'], [0.0])

        return loadCurveDict

    def loadCurveFromDict(self, lcid, lcDict, sfo=1.0, sfa=1.0, offo=0.0, offa=0.0, dattyp=0, lcint=0, shearHeight=1.0,
                          transform=None):
        '''shape = 'vel' accelerates in endTime/100'''
        # offa offsets in time, to the right

        sidr = 0
        lcint = 0

        line_block = ['*DEFINE_CURVE_TITLE\n']
        line_block.append('Load curve {}\n'.format(lcid))
        line_block.append(self.format_comment_line(['lcid', 'sidr', 'sfa', 'sfo', 'offa', 'offo', 'dattyp', 'lcint']))
        line_block.append(self.format_key_line([lcid, sidr, sfa, sfo, offa, offo, dattyp, lcint]))
        line_block.append(self.format_comment_line(['a1', 'o1']))
        if transform == None:
            for abscissa, ordinate in zip(lcDict['time'], lcDict['velAmpl']):
                line_block.append(self.format_key_line([abscissa, ordinate]))
        if transform == 'shearAngle':
            for abscissa, ordinate in zip(lcDict['time'], lcDict['disp']):
                line_block.append(self.format_key_line([abscissa, np.arctan(ordinate / shearHeight)]))
        self.submit_block(line_block)

    def loadCurve(self, lcid, startTime=0.0, endTime=1.0, shape='linear', sfo=1.0, offo=0.0, tanRatioEnd=0.0):
        '''shape = 'vel' accelerates in endTime/100'''
        # offa offsets in time, to the right

        sidr = 0
        sfa = 1.0
        offa = 0.0
        dattyp = 0
        lcint = 0
        endAmp = 1.0

        line_block = ['*DEFINE_CURVE_TITLE\n']
        line_block.append('Load curve {}\n'.format(lcid))
        line_block.append(self.format_comment_line(['lcid', 'sidr', 'sfa', 'sfo', 'offa', 'offo', 'dattyp', 'lcint']))
        line_block.append(self.format_key_line([lcid, sidr, sfa, sfo, offa, offo, dattyp, lcint]))
        line_block.append(self.format_comment_line(['a1', 'o1']))
        line_block.append(self.format_key_line([0.0, 0.0]))
        if startTime != 0.0:
            line_block.append(self.format_key_line([0.0, 0.0]))
            line_block.append(self.format_key_line([startTime, 0.0]))
        if shape == 'vel':
            line_block.append(self.format_key_line([endTime / 100, endAmp]))
            line_block.append(self.format_key_line([endTime, endAmp]))
            line_block.append(self.format_key_line([endTime + 100, endAmp]))
        elif shape == 'tan':
            abscissaArray = np.linspace(startTime, endTime, 1000)[1:-1]
            ordinateValue = np.arctan(np.linspace(0, tanRatioEnd, 1000))[1:-1]
            for absc, ordin in zip(abscissaArray, ordinateValue):
                line_block.append(self.format_key_line([absc, ordin]))
            line_block.append(self.format_key_line([endTime, np.arctan(tanRatioEnd)]))
            line_block.append(self.format_key_line([endTime + 100, np.arctan(tanRatioEnd)]))
        elif shape == 'linear':
            line_block.append(self.format_key_line([endTime, endAmp]))
            line_block.append(self.format_key_line([endTime + 100, endAmp]))
        elif shape == 'flat':
            endAmp = 0.0
            line_block.append(self.format_key_line([endTime, endAmp]))
            line_block.append(self.format_key_line([endTime + 100, endAmp]))

        self.submit_block(line_block)

    ###########################################
    # BCs
    ################################################
    def constrained_multiple_global(self, id, constr_list, direction):
        line_block = ['*CONSTRAINED_MULTIPLE_GLOBAL\n']
        line_block.append(self.format_comment_line(['id']))
        line_block.append(self.format_key_line_short([id]))
        for const_dict in constr_list:
            line_block.append(self.format_key_line_short([const_dict['nmp']]))
            # line_block.append(self.format_comment_line([     nid       dir      coef']))
            for nid, coeff in zip(const_dict['nid_list'], const_dict['coeff_list']):
                line_block.append(self.format_key_line_short([nid, direction, float(coeff)]))
        self.submit_block(line_block)

    def constrained_linear_local(self, lcid, nid_list, coeff_list, direction, cid=0):
        line_block = ['*CONSTRAINED_LINEAR_LOCAL\n']
        line_block.append(self.format_comment_line(['lcid']))
        line_block.append(self.format_key_line([lcid]))
        line_block.append(self.format_comment_line(['nid', 'dir', 'cid', 'coef']))
        for nid, coeff in zip(nid_list, coeff_list):
            line_block.append(self.format_key_line([nid, direction, cid, float(coeff)]))
        self.submit_block(line_block)

    def constrained_linear_global(self, lcid, nid_list, coef_list, direction):
        line_block = ['*CONSTRAINED_LINEAR_GLOBAL\n']
        line_block.append(self.format_comment_line(['lcid']))
        line_block.append(self.format_key_line([lcid]))
        line_block.append(self.format_comment_line(['nid', 'dir', 'coef']))
        for nid, coef in zip(nid_list, coef_list):
            line_block.append(self.format_key_line([nid, direction, float(coef)]))
        self.submit_block(line_block)

    def rigid_wall_planar(self, rwid, norm_vect, fric=0.01, offset=0.0):
        nsid = 0
        nsidex = 0
        boxid = 0
        rwsf = 1.0
        wvel = 0.0

        vectTail = norm_vect[0]
        vectHead = norm_vect[1]
        line_block = ['*RIGIDWALL_PLANAR_ID\n']
        line_block.append(self.format_key_line([rwid, 'RW planar {}'.format(rwid)]))
        line_block.append(self.format_comment_line(['nsid', 'nsidex', 'boxid', 'offset', 'birth', 'death', 'rwksf']))
        line_block.append(self.format_key_line([nsid, nsidex, boxid, offset, 0.0, 1.0e20, rwsf]))
        line_block.append(self.format_comment_line(['xt', 'yt', 'zt', 'xh', 'yh', 'zh', 'fric', 'wvel']))
        line_block.append(self.format_key_line(
            [vectTail[0], vectTail[1], vectTail[2], vectHead[0], vectHead[1], vectHead[2], fric, wvel]))
        self.submit_block(line_block)

    def rigid_wall_geometric_flatt_motion(self, rwid, lcid, norm_vect, vect_L_head, vel_vect, fric=0.0, opt=1):
        nsid = 0
        nsidex = 0
        boxid = 0
        lenl = 0.0
        lenm = 0.0
        vect_tail = norm_vect[0]
        vect_head = norm_vect[1]
        line_block = ['*RIGIDWALL_GEOMETRIC_FLAT_MOTION_ID\n']
        line_block.append(self.format_key_line([rwid, 'RW Flat Motion {}'.format(rwid)]))
        line_block.append(self.format_comment_line(['nsid', 'nsidex', 'boxid', 'birth', 'death']))
        line_block.append(self.format_key_line([nsid, nsidex, boxid, 0.0, 1.0e20]))
        line_block.append(self.format_comment_line(['xt', 'yt', 'zt', 'xh', 'yh', 'zh', 'fric']))
        line_block.append(
            self.format_key_line([vect_tail[0], vect_tail[1], vect_tail[2], vect_head[0], vect_head[1], vect_head[2], fric]))
        line_block.append(self.format_comment_line(['xhev', 'yhev', 'zhev', 'lenl', 'lenm']))
        line_block.append(self.format_key_line([vect_L_head[0], vect_L_head[1], vect_L_head[2], lenl, lenm]))
        line_block.append(self.format_comment_line(['lcid', 'opt', 'vx', 'vy', 'vz']))
        line_block.append(self.format_key_line([lcid, opt, vel_vect[0], vel_vect[1], vel_vect[2]]))
        self.submit_block(line_block)

    def boundary_spc_set(self, bspcid, nsid, dofx=1, dofy=1, dofz=1, dofrx=0, dofry=0, dofrz=0, cid=0):
        line_block = ['*BOUNDARY_SPC_SET_ID\n']
        line_block.append(self.format_key_line([bspcid, 'Boundary SPC {}'.format(bspcid)]))
        line_block.append(self.format_comment_line(['nsid', 'cid', 'dofx', 'dofy', 'dofz', 'dofrx', 'dofry', 'dofrz']))
        line_block.append(self.format_key_line([nsid, cid, dofx, dofy, dofz, dofrx, dofry, dofrz]))
        self.submit_block(line_block)

    def boundary_spc_node(self, bspcid, nid, dofx=1, dofy=1, dofz=1, dofrx=0, dofry=0, dofrz=0, cid=0):
        line_block = ['*BOUNDARY_SPC_NODE_ID\n']
        line_block.append(self.format_key_line([bspcid, 'Boundary SPC {}'.format(bspcid)]))
        line_block.append(self.format_comment_line(['nsid', 'cid', 'dofx', 'dofy', 'dofz', 'dofrx', 'dofry', 'dofrz']))
        line_block.append(self.format_key_line([nid, cid, dofx, dofy, dofz, dofrx, dofry, dofrz]))
        self.submit_block(line_block)

    def boundary_prescribed_motion_set(self, bmsid, nsid, dof, vad, lcid, sf=1.0, death=1.0E28, birth=0.0, offset1=0.0,
                                       offset2=0.0):
        # dof: x, y, x= 1, 2, 3
        vid = 0
        mrb = 0
        node1 = 0
        node2 = 0
        line_block = ['*BOUNDARY_PRESCRIBED_MOTION_SET_ID\n']
        line_block.append(self.format_key_line([bmsid, 'Boundary disp {}'.format(bmsid)]))
        line_block.append(self.format_comment_line(['nsid', 'dof', 'vad', 'lcid', 'sf', 'vid', 'death', 'birth']))
        line_block.append(self.format_key_line([nsid, dof, vad, lcid, sf, vid, death, birth]))
        if abs(dof) in [9, 10, 11]:
            line_block.append(self.format_comment_line(['offset1', 'offset2', 'mrb', 'node1', 'node2']))
            line_block.append(self.format_key_line([offset1, offset2, mrb, node1, node2]))

        self.submit_block(line_block)

    def boundary_prescribed_motion_node(self, bmsid, nid, dof, vad, lcid, sf=1.0, death=1.0E28, birth=0.0, offset1=0.0,
                                        offset2=0.0):
        # dof: x, y, x= 1, 2, 3
        vid = 0
        mrb = 0
        node1 = 0
        node2 = 0
        line_block = ['*BOUNDARY_PRESCRIBED_MOTION_NODE_ID\n']
        line_block.append(self.format_key_line([bmsid, 'Boundary disp {}'.format(bmsid)]))
        line_block.append(self.format_comment_line(['nsid', 'dof', 'vad', 'lcid', 'sf', 'vid', 'death', 'birth']))
        line_block.append(self.format_key_line([nid, dof, vad, lcid, sf, vid, death, birth]))
        if abs(dof) in [9, 10, 11]:
            line_block.append(self.format_comment_line(['nsid', 'dof', 'vad', 'lcid', 'sf', 'vid', 'death', 'birth']))
            line_block.append(self.format_comment_line(['offset1', 'offset2', 'mrb', 'node1', 'node2']))
            line_block.append(self.format_key_line([offset1, offset2, mrb, node1, node2]))

        self.submit_block(line_block)

    ###############################################
    # Contact
    ##############################################\
    def contact_automatic_general_interior_id(self, cid, ssid=0, sstyp=0, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=0):
        line_block = ['*CONTACT_AUTOMATIC_GENERAL_INTERIOR_ID\n']
        line_block.append(self.format_key_line([cid, 'Automatic contact full']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, 0, sstyp, 0, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))

        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
        line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
        line_block.append(self.format_key_line([0.0000000000e+000, 1, 1, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
        line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
        line_block.append(
            self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))

        self.submit_block(line_block)

    def contact_automatic_general_id(self, cid, ssid=0, sstyp=0, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=0):
        line_block = ['*CONTACT_AUTOMATIC_GENERAL_ID\n']
        line_block.append(self.format_key_line([cid, 'Automatic contact full']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, 0, sstyp, 0, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))

        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
        line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
        line_block.append(self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
        line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
        line_block.append(
            self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))

        self.submit_block(line_block)

    def contact_single_surface_id(self, cid, ssid=0, sstyp=0, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=1, depth=1, igap=0):
        line_block = ['*CONTACT_SINGLE_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'SingleSurface']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, 0, sstyp, 0, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))

        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 3.0, depth, 0, 1.0]))
        line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
        line_block.append(
            self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
        line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
        line_block.append(
            self.format_key_line([igap, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))

        self.submit_block(line_block)

    def contact_automatic_single_surface_id(self, cid, ssid=0, sstyp=0, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=1, depth=1,
                                            igap=1, snlog=1):
        line_block = ['*CONTACT_AUTOMATIC_SINGLE_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'Automatic contact full']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, 0, sstyp, 0, 0, 0, 1, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))

        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, depth, 0, 1.0]))
        line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
        line_block.append(
            self.format_key_line([0.0000000000e+000, 0, 0, snlog, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
        line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
        line_block.append(
            self.format_key_line([igap, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
        self.submit_block(line_block)

    def contact_automatic_single_surface_mortar_id(self, cid, ssid=0, sstyp=0, fs=-1.0, fd=-1.0, dc=-1.0, soft=0, ignore=0,
                                                   depth=3, igap=1, snlog=0, vdc=0):
        line_block = ['*CONTACT_AUTOMATIC_SINGLE_SURFACE_MORTAR_ID\n']
        line_block.append(self.format_key_line([cid, 'Automatic contact full']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, 0, sstyp, 0, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, vdc, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))

        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, depth, 0, 1.0]))
        line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
        line_block.append(
            self.format_key_line([0.0000000000e+000, 0, 0, snlog, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
        line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
        line_block.append(
            self.format_key_line([igap, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
        self.submit_block(line_block)

    def contact_automatic_nodes_to_surface(self, cid, ssid, msid, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=1, depth=1):
        line_block = ['*CONTACT_AUTOMATIC_NODES_TO_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'NodesToSurface contact']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, 4, 3, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))
        if ignore != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
            line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
            line_block.append(
                self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
            line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
            line_block.append(
                self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
        elif soft != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, depth, 0, 1.0]))
        self.submit_block(line_block)

    def contact_tied_nodes_to_surface_id(self, cid, ssid, msid, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=1):
        line_block = ['*CONTACT_TIED_NODES_TO_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'Tied contact']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, 4, 0, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))
        if ignore != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
            line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
            line_block.append(
                self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
            line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
            line_block.append(
                self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
        elif soft != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
        self.submit_block(line_block)

    def contact_tied_shell_edge_to_surface_offset_id(self, cid, ssid, msid, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=0):
        line_block = ['*CONTACT_TIED_SHELL_EDGE_TO_SURFACE_OFFSET_ID\n']
        line_block.append(self.format_key_line([cid, 'Tied contact']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, 4, 3, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))
        if ignore != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
            line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
            line_block.append(
                self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
            line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
            line_block.append(
                self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
        elif soft != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
        self.submit_block(line_block)

    def contact_tied_shell_edge_to_surface_id(self, cid, ssid, msid, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=0):
        line_block = ['*CONTACT_TIED_SHELL_EDGE_TO_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'Tied contact']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, 2, 3, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))
        if ignore != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
            line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
            line_block.append(
                self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
            line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
            line_block.append(
                self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
        elif soft != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
        self.submit_block(line_block)

    def contact_constraint_nodes_to_surface(self, cid, ssid, msid, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=0):
        line_block = ['*CONTACT_CONSTRAINT_NODES_TO_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'Tied contact']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, 4, 0, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))
        line_block.append(self.format_comment_line(['kpf']))
        line_block.append(self.format_key_line([1.0]))
        if ignore != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
            line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
            line_block.append(
                self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
            line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
            line_block.append(
                self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
        elif soft != 0:
            line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
            line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
        self.submit_block(line_block)

    def contact_automatic_surface_to_surface_tie_break(self, cid, ssid, msid, fs=0.0, fd=0.0, dc=0.0, soft=1, ignore=0,
                                                       option=2, tiedid=0, shledg=0):
        line_block = ['*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE_TIEBREAK_ID\n']
        line_block.append(self.format_key_line([cid, 'Tiebreak contact']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, 2, 3, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))
        line_block.append(self.format_comment_line(['option', 'nfls', 'sfls', 'param', 'eraten', 'erates', 'ct2cn', 'cn']))
        line_block.append(self.format_key_line([option, 10E8, 10E8, 0, 0, 0, 0, 0]))  # alt option 4

        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))
        if 1 == 1:
            igap = 0
            line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
            line_block.append(
                self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
            line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
            line_block.append(self.format_key_line(
                [igap, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
            line_block.append(self.format_comment_line(['q2tri', 'dtpchk', 'sfnbr', 'fnlscl', 'dnlscl', 'tcso', 'tiedid', 'shledg']))
            line_block.append(self.format_key_line([0, 0.0, 0.0, 0., 0., 0, tiedid, shledg]))
        self.submit_block(line_block)

    def contact_tiebreak_nodes_to_surface(self, cid, ssid, msid, fs=0.0, fd=0.0, dc=0.0, soft=1, ignore=0, option=2,
                                          tiedid=0, shledg=0):
        line_block = ['*CONTACT_TIEBREAK_NODES_TO_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'Tiebreak contact']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, 4, 3, 0, 0, 0, 1]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))
        line_block.append(self.format_comment_line(['nflf', 'sflf', 'nen', 'mes']))
        line_block.append(self.format_key_line([10E7, 10E7, 2, 2]))  # alt option 4
        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 2.0, 2.0, 0, 1.0]))

        igap = 0
        line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
        line_block.append(
            self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
        line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
        line_block.append(self.format_key_line(
            [igap, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))
        line_block.append(self.format_comment_line(['q2tri', 'dtpchk', 'sfnbr', 'fnlscl', 'dnlscl', 'tcso', 'tiedid', 'shledg']))
        line_block.append(self.format_key_line([0, 0.0, 0.0, 0., 0., 0, tiedid, shledg]))
        self.submit_block(line_block)

    def contact_force_transducer(self, cid, ssid):
        line_block = ['*CONTACT_FORCE_TRANSDUCER_PENALTY_ID\n']
        line_block.append(self.format_comment_line(['cid']))
        line_block.append(self.format_key_line([cid, 'Force Transducer']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, 0, 3, 0, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))

        self.submit_block(line_block)

    def contact_surface_to_surface_id(self, cid, ssid=0, msid=0, sstyp=4, mstyp=3, fs=0.0, fd=0.0, dc=0.0, soft=0, ignore=1,
                                      depth=1):
        line_block = ['*CONTACT_SURFACE_TO_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'SingleSurface']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, sstyp, mstyp, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))

        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 3.0, depth, 0, 1.0]))
        line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
        line_block.append(
            self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
        line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
        line_block.append(
            self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))

        self.submit_block(line_block)

    def contact_airbag_single_surface_id(self, cid, ssid=0, msid=0, sstyp=4, mstyp=3, fs=0.0, fd=0.0, dc=0.0, soft=0,
                                         ignore=1, depth=1):
        line_block = ['*CONTACT_AIRBAG_SINGLE_SURFACE_ID\n']
        line_block.append(self.format_key_line([cid, 'SingleSurface']))
        line_block.append(self.format_comment_line(['ssid', 'msid', 'sstyp', 'mstyp', 'sboxid', 'mboxid', 'spr', 'mpr']))
        line_block.append(self.format_key_line([ssid, msid, sstyp, mstyp, 0, 0, 0, 0]))
        line_block.append(self.format_comment_line(['fs', 'fd', 'dc', 'vc', 'vdc', 'penchk', 'bt', 'dt']))
        line_block.append(self.format_key_line([fs, fd, dc, 0.0, 0.0, 0, 0.0, 1.00000E20]))
        line_block.append(self.format_comment_line(['sfs', 'sfm', 'sst', 'mst', 'sfst', 'sfmt', 'fsf', 'vsf']))
        line_block.append(self.format_key_line([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]))

        line_block.append(self.format_comment_line(['soft', 'sofscl', 'lcidab', 'maxpar', 'sbopt', 'depth', 'bsort', 'frcfrq']))
        line_block.append(self.format_key_line([soft, 0.1, 0, 1.025, 3.0, depth, 0, 1.0]))
        line_block.append(self.format_comment_line(['penmax', 'thkopt', 'shlthk', 'snlog', 'isym', 'i2d3d', 'sldthk', 'sldstf']))
        line_block.append(
            self.format_key_line([0.0000000000e+000, 0, 0, 0, 0, 0, 0.0000000000e+000, 0.0000000000e+000]))
        line_block.append(self.format_comment_line(['igap', 'ignore        dprfac/mpar1        dtstif/mpar2', 'unused', 'unused', 'flangl', 'cid_rcf']))
        line_block.append(
            self.format_key_line([1, ignore, 0.0000000000e+000, 0.0000000000e+000, ' ', ' ', 0.0000000000e+000, 0]))

        self.submit_block(line_block)

    def pertubation_node(self, scl, nsid=0, cmp=7, xwl=1.0, ywl=1.0, zwl=1.0):
        type = 1
        ampl = 1
        xoff = 0.0
        yoff = 0.0
        zoff = 0.0
        line_block = ['*PERTURBATION_NODE\n']
        line_block.append(self.format_comment_line(['type', 'nsid', 'scl', 'cmp', 'icoord', 'cid']))
        line_block.append(self.format_key_line([type, nsid, scl, cmp, 0, 0]))
        line_block.append(self.format_comment_line(['ampl', 'xwl', 'xoff', 'ywl', 'yoff', 'zwl', 'zoff']))
        line_block.append(self.format_key_line([ampl, xwl, xoff, ywl, yoff, zwl, zoff]))
        self.submit_block(line_block)

    def pertubation_shell(self, scl, nsid=0, cmp=7, xwl=1.0, ywl=1.0, zwl=1.0):
        type = 1
        ampl = 1
        xoff = 0.0
        yoff = 0.0
        zoff = 0.0
        line_block = ['*PERTURBATION_SHELL_THICKNESS\n']
        line_block.append(self.format_comment_line(['type', 'nsid', 'scl', 'cmp', 'icoord', 'cid']))
        line_block.append(self.format_key_line([type, nsid, scl, cmp, 0, 0]))
        line_block.append(self.format_comment_line(['ampl', 'xwl', 'xoff', 'ywl', 'yoff', 'zwl', 'zoff']))
        line_block.append(self.format_key_line([ampl, xwl, xoff, ywl, yoff, zwl, zoff]))
        self.submit_block(line_block)

    def airbag_simple_pressure_volume(self, abid, sid):
        sidtyp = 1
        rbid = 0
        vsca = 1.0
        psca = 1.0
        vini = 0.0
        mwd = 0.0
        spsf = 0.0
        cn = 0.1
        beta = 1.0
        lcid = 0
        lciddr = 0
        line_block = ['*AIRBAG_SIMPLE_PRESSURE_VOLUME_ID\n']
        line_block.append(self.format_comment_line(['abid']))
        line_block.append(self.format_key_line([abid]))
        line_block.append(self.format_comment_line(['sid', 'sidtyp', 'rbid', 'vsca', 'psca', 'vini', 'mwd', 'spsf']))
        line_block.append(self.format_key_line([sid, sidtyp, rbid, vsca, psca, vini, mwd, spsf]))
        line_block.append(self.format_comment_line(['cn', 'beta', 'lcid', 'lciddr']))
        line_block.append(self.format_key_line([cn, beta, lcid, lciddr]))
        self.submit_block(line_block)

    def airbag_adiabatic_gas_model(self, abid, sid):
        sidtyp = 1
        rbid = 0
        vsca = 1.0
        psca = 1.0
        vini = 0.0
        mwd = 0.0
        spsf = 0.0
        psf = 1.0
        lcid = 0
        gamma = 1.4
        po = 0.0
        pa = 0.1
        ro = 1.3e-12
        line_block = ['*AIRBAG_ADIABATIC_GAS_MODEL_ID\n']
        line_block.append(self.format_comment_line(['abid']))
        line_block.append(self.format_key_line([abid, 'Adiabatic airbag']))
        line_block.append(self.format_comment_line(['sid', 'sidtyp', 'rbid', 'vsca', 'psca', 'vini', 'mwd', 'spsf']))
        line_block.append(self.format_key_line([sid, sidtyp, rbid, vsca, psca, vini, mwd, spsf]))
        line_block.append(self.format_comment_line(['psf', 'lcid', 'gamma', 'po', 'pe', 'ro']))
        line_block.append(self.format_key_line([psf, lcid, gamma, po, pa, ro]))
        self.submit_block(line_block)
