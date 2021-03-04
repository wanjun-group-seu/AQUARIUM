# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import os

import pysrc.body.logger

__doc__ = '''
'''
__author__ = 'zerodel'

_logger = pysrc.body.logger.default_logger("opt_check")

_description_template_default = """


[{working}]
# this section [{working}] contains following options:
#######  essential arguments

{must_have}

#######  only one of these option pairs are needed

{equivalent}

#######  optional arguments

{optional}

#######  only one of these following options are optional

{equivalent_optional}


#######  these options are not allowed

{banned}


        """


class OptionChecker(object):
    def __init__(self, opt_dict=None, name="", logger=_logger):
        self.name = name
        self.dict_opt = opt_dict
        self.opts = []
        self.opts_equivalent = []
        self.opts_forbidden = []
        self.opts_optional = []
        self.opts_optional_equivalent = []

        self.ad_hoc_conditions = []

        self.check_func = {}
        self.exceptions = {}
        self.descriptions = {}
        self._logger = logger

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None and exc_val is None and exc_tb is None:
            self.check()

    def _self_description(self):
        # todo: this should have some modification to be put into standard error .
        print(str(self))

    def check(self, opt_dict=None):
        if opt_dict:
            self.dict_opt = opt_dict

        try:
            for opt_str in self.opts:
                self._check_opt(opt_str)
            for opt_str in self.opts_forbidden:
                if opt_str in self.dict_opt:
                    raise KeyError("Error: This argument is forbidden: {}".format(opt_str))
            for opt_str in self.opts_optional:
                if opt_str in self.dict_opt:
                    self._check_opt(opt_str)
            for opt_group in self.opts_equivalent:
                self._check_equivalent_arg_list(opt_group)

            for opt_group in self.opts_optional_equivalent:
                self._check_optional_equivalent_arg_list(opt_group)

            for cond_and_desc in self.ad_hoc_conditions:
                cond_and_desc[0](self.dict_opt)

        except Exception as e:
            self._self_description()
            print(e)
            raise e

    def _lazy_add(self, opt_str, check_fun, exception_on_error, desc=""):
        self.opts_optional.append(opt_str)
        self.check_func[opt_str] = check_fun
        self.exceptions[opt_str] = exception_on_error
        self.descriptions[opt_str] = desc

    def _add_strict(self, opt_str, check_fun, exception_on_error, desc=""):
        self.opts.append(opt_str)
        self.check_func[opt_str] = check_fun
        self.exceptions[opt_str] = exception_on_error
        self.descriptions[opt_str] = desc

    def must_have(self, opt_str, check_fun, exception_on_error, des_str=""):
        self._add_strict(opt_str, check_fun, exception_on_error, des_str)

    def may_need(self, opt_str, check_fun, exception_on_error, des_str=""):
        self._lazy_add(opt_str, check_fun, exception_on_error, des_str)

    def one_and_only_one(self, equivalent_args, check_fun, exception_on_error, des_str=""):
        if isinstance(equivalent_args, list):
            self.opts_equivalent.append(equivalent_args)
            self.descriptions[tuple(equivalent_args)] = des_str
            for opt in equivalent_args:
                self.check_func[opt] = check_fun
                self.exceptions[opt] = exception_on_error
                self.descriptions[opt] = des_str
        else:
            raise KeyError("Error@OptionChecker: need a list of equivalent args, which should have at least 2 args")

    def at_most_one(self, equivalent_args, check_fun, exception_on_error, des_str=""):
        if isinstance(equivalent_args, list):
            self.opts_optional_equivalent.append(equivalent_args)
            self.descriptions[tuple(equivalent_args)] = des_str
            for opt in equivalent_args:
                self.check_func[opt] = check_fun
                self.exceptions[opt] = exception_on_error
                self.descriptions[opt] = des_str
        else:
            raise KeyError("Error@OptionChecker: need a list of equivalent args, which should have at least 2 args")

    def _check_optional_equivalent_arg_list(self, optional_equivalent_args):
        args_found = [x for x in optional_equivalent_args if x in self.dict_opt]
        num_args_found = len(args_found)

        self._logger.debug("checking equivalent opts :{}".format(" ".join(optional_equivalent_args)))

        self._logger.debug("[ {} ] found in equivalent opts : [ {} ]".format(" ".join(args_found),
                                                                             " ".join(optional_equivalent_args)))

        if num_args_found > 1:
            raise KeyError("ERROR@OptionChecker: redundant argument found, {}".format("/".join(args_found)))
        elif num_args_found == 1:
            that_only_args = args_found[0]
            self._check_the_only_one_in_equivalent_args(that_only_args)
        else:  # no arguments found at command line .
            pass

    def _check_equivalent_arg_list(self, equivalent_args):
        args_in_opts = [x for x in equivalent_args if x in self.dict_opt]
        num_args = len(args_in_opts)

        self._logger.debug("checking equivalent opts : [ %s ]" % " ".join(equivalent_args))

        self._logger.debug("[ {opt} ] found in equivalent opts : [ {opt_group} ]".format(opt=" ".join(args_in_opts),
                                                                                         opt_group=" ".join(
                                                                                             equivalent_args)))

        if num_args == 1:
            that_option = args_in_opts[0]

            self._logger.debug("{} is the only one opt in [ {} ]".format(that_option,
                                                                         " ".join(equivalent_args)))

            self._check_the_only_one_in_equivalent_args(that_option)

        elif num_args == 0:
            raise KeyError("Error@OptionChecker: lack of essential option : {} ".format("/".join(equivalent_args)))
        else:
            raise KeyError("Error@OptionChecker: Redundant parameter : {}".format("/".join(args_in_opts)))

    def _check_the_only_one_in_equivalent_args(self, that_option):

        that_value_of_option = self.dict_opt[that_option]
        check_fun = self.check_func[that_option]
        exception_on_error = self.exceptions[that_option]
        if not check_fun(that_value_of_option):
            raise exception_on_error

    def _check_opt(self, opt_str):

        self._logger.debug("checking %s" % opt_str)

        if not isinstance(opt_str, str):
            raise TypeError("Error@OptionChecker: only accept string input : {}".format(str(opt_str)))

        if opt_str not in self.dict_opt:
            raise KeyError("Error@OptionCheck: no such option {}, unable to perform check".format(opt_str))

        val_this_opt = self.dict_opt[opt_str]
        self._logger.debug("checking {key} : value of {key} is {val}".format(key=opt_str,
                                                                             val=str(val_this_opt)))

        if not self.check_func[opt_str](val_this_opt):
            _logger.error("Error@OptionCheck: fail on key : {key} value: {val}".format(key=opt_str, val=val_this_opt))
            raise self.exceptions[opt_str]

    def forbid_these_args(self, *args):
        self.opts_forbidden.extend(args)

    def __str__(self):

        def _add_sharp_head(str_description):
            return "\n".join(["# %s" % x.strip() for x in str_description.strip().split("\n")])

        def _inner_make_pair_opt_desc(opts, desc_dict=self.descriptions):
            out = []
            for opt in opts:
                str_description = _add_sharp_head(desc_dict.get(opt, ""))
                if self.dict_opt and opt in self.dict_opt:
                    out.append(
                        "{opt} = {val}\n{desc}".format(opt=opt, desc=str_description, val=self.dict_opt.get(opt, "")))
                else:
                    out.append("{opt}\n{desc}".format(opt=opt, desc=str_description))
            return out

        def _display_option_group(option_groups, description_dict):
            return "\n\n".join(
                ["%s\n%s" % (
                    '/'.join(grp), _add_sharp_head(description_dict[tuple(grp)])) for
                 grp
                 in option_groups])

        must_have = "\n\n".join(_inner_make_pair_opt_desc(self.opts, self.descriptions))

        optional = "\n\n".join(_inner_make_pair_opt_desc(self.opts_optional, self.descriptions))

        banned = "\n\n".join(["# # %s" % x.strip() for x in self.opts_forbidden]) if isinstance(self.opts_forbidden,
                                                                                                list) else "# %s" % self.opts_forbidden

        equivalent = _display_option_group(self.opts_equivalent, self.descriptions)

        optional_equivalent = _display_option_group(self.opts_optional_equivalent, self.descriptions)

        return _description_template_default.format(working=self.name,
                                                    must_have=must_have,
                                                    equivalent=equivalent,
                                                    optional=optional,
                                                    equivalent_optional=optional_equivalent,
                                                    banned=banned)

    def custom_condition(self, custom_check_fun, desc=""):
        self.ad_hoc_conditions.append((custom_check_fun, desc))


def this_abs_path_should_be_a_dir(path_dir_you_want):
    if os.path.exists(path_dir_you_want):
        _logger.debug("already has a file/dir in %s" % path_dir_you_want)
        return os.path.isdir(path_dir_you_want)
    else:
        _logger.debug("no directory at %s, make a new one " % path_dir_you_want)
        os.makedirs(path_dir_you_want, exist_ok=True)


def all_files_in_this_string_should_exist(file_list_str):
    if file_list_str:
        file_parts = file_list_str.strip().split()
        for single_file in file_parts:
            _logger.debug("checking file %s" % single_file)
            if not os.path.isfile(single_file):
                _logger.error("error occurs when checking %s" % single_file)
                return False
        else:
            return True
    else:
        _logger.warning("suppose to be a list of file , nothing found")
        return False
