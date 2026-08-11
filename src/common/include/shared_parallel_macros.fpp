#:set NVIDIA_COMPILER_ID="NVHPC"
#:set PGI_COMPILER_ID="PGI"
#:set INTEL_COMPILER_ID="Intel"
#:set CCE_COMPILER_ID="Cray"
#:set AMD_COMPILER_ID="LLVMFlang"

#:set USING_NVHPC = (MFC_COMPILER == NVIDIA_COMPILER_ID or MFC_COMPILER == PGI_COMPILER_ID)
#:set USING_CCE = (MFC_COMPILER == CCE_COMPILER_ID)
#:set USING_AMD = (MFC_COMPILER == AMD_COMPILER_ID)

#:def ASSERT_LIST(data, datatype)
    #:assert data is not None
    #:assert isinstance(data, list)
    #:assert len(data) != 0
    #:assert all(isinstance(element, datatype) for element in data)
#:enddef

#:def GEN_CLAUSE(clause_name, clause_str)
    #:set clause_regex = re.compile(',(?![^(]*\\))')
    #:assert isinstance(clause_name, str)
    #:if clause_str is not None
        #:set count = 0
        #:assert isinstance(clause_str, str)
        #:assert clause_str[0] == '[' and clause_str[-1] == ']'
        #:for c in clause_str
            #:if c == '('
                #:set count = count + 1
            #:elif c == ')'
                #:set count = count - 1
            #:endif
            #:if c == ',' and count > 1
                #:stop 'Nested parentheses with comma inside is not supported. Incorrect clause: {}'.format(clause_str)
            #:elif count < 0
                #:stop 'Missing parentheses. Incorrect clause: {}'.format(clause_str)
            #:endif
        #:endfor
        #:set clause_str = re.sub(clause_regex, ';', clause_str)
        #:set clause_list = [x.strip() for x in clause_str.strip('[]').split(';')]
        $:ASSERT_LIST(clause_list, str)
        #:set clause_str = clause_name + ', '.join(clause_list) + ') '
    #:else
        #:set clause_str = ''
    #:endif
    $:clause_str
#:enddef

#:def GEN_PARENTHESES_CLAUSE(clause_name, clause_str)
    #:assert isinstance(clause_name, str)
    #:if clause_str is not None
        #:assert isinstance(clause_str, str)
        #:set clause = clause_name + '('
        #:set clause_str = GEN_CLAUSE(clause, clause_str)
    #:else
        #:set clause_str = ''
    #:endif
    $:clause_str
#:enddef

#:def GEN_PRIVATE_STR(private, initialized_values)
    #:assert isinstance(initialized_values, bool)
    #:if initialized_values == True
        #:set private_val = GEN_PARENTHESES_CLAUSE('firstprivate', private)
    #:else
        #:set private_val = GEN_PARENTHESES_CLAUSE('private', private)
    #:endif
    $:private_val
#:enddef

#:def GEN_REDUCTION_STR(reduction, reductionOp)
    #:if reduction is not None and reductionOp is not None
        #:assert isinstance(reduction, str)
        #:assert isinstance(reductionOp, str)
        #:assert reduction[0] == '[' and reduction[-1] == ']'
        #:assert reductionOp[0] == '[' and reductionOp[-1] == ']'
        #:set reduction = reduction.replace(' ', '')
        #:set reduction = reduction[1:-1]
        #:set reduction_list = reduction.split('],')
        #:set reduction_list = [str + ']' for str in reduction_list]
        #:assert all(str[0] == '[' and str[-1] == ']' for str in reduction_list)

        #:set reductionOp_list = [x.strip() for x in reductionOp.strip('[]').split(',')]
        $:ASSERT_LIST(reduction_list, str)
        $:ASSERT_LIST(reductionOp_list, str)
        #:assert len(reduction_list) == len(reductionOp_list)
        #:set reduction_val = ''
        #:for i in range(len(reduction_list))
            #:set temp_clause = GEN_PARENTHESES_CLAUSE('reduction', reduction_list[i]).strip('\n')
            #:set ind = temp_clause.find('reduction(') + len('reduction(')
            #:set reduction_val = reduction_val.strip('\n') + temp_clause[:ind] + reductionOp_list[i] + ':' + temp_clause[ind:]
        #:endfor
    #:elif reduction is not None or reductionOp is not None
        #:stop 'Cannot set the reduction list or reduction operation without setting the other'
    #:else
        #:set reduction_val = ''
    #:endif
    $:reduction_val
#:enddef

#:def GEN_COLLAPSE_STR(collapse)
    #:if collapse is not None
        #:set collapse = int(collapse)
        #:assert isinstance(collapse, int)
        #:assert collapse > 1
        #:set collapse_val = 'collapse(' + str(collapse) + ') '
    #:else
        #:set collapse_val = ''
    #:endif
    $:collapse_val
#:enddef

#:def GEN_EXTRA_ARGS_STR(extraArgs)
    #:if extraArgs is not None
        #:assert isinstance(extraArgs, str)
        #:set extraArgs_val = extraArgs
    #:else
        #:set extraArgs_val = ''
    #:endif
    $:extraArgs_val
#:enddef

#:def FOLD_DIRECTIVE(directive, sentinel, width=200)
    #! Fold a long GPU directive across free-form continuation lines so it stays under
    #! nvfortran's ~1000-char source-line limit, breaking at whole-clause boundaries
    #! (clause(args) groups and bare keywords) and repeating the sentinel (e.g. '!$acc&') on
    #! each continuation -- which fypp's --no-folding cannot do (its generic folder omits the
    #! sentinel). A single clause wider than `width` is split at its top-level commas too.
    #:set _clauses = re.findall(r'\w+\([^)]*\)|\S+', directive)
    #:set _toks = []
    #:for _cl in _clauses
        #:if len(_cl) > width and ',' in _cl
            #:set _buf = ''
            #:set _depth = 0
            #:for _ch in _cl
                #:if _ch == '('
                    #:set _depth = _depth + 1
                #:elif _ch == ')'
                    #:set _depth = _depth - 1
                #:endif
                #:set _buf = _buf + _ch
                #:if _ch == ',' and _depth == 1
                    #:set _ = _toks.append(_buf.strip())
                    #:set _buf = ''
                #:endif
            #:endfor
            #:if _buf.strip() != ''
                #:set _ = _toks.append(_buf.strip())
            #:endif
        #:else
            #:set _ = _toks.append(_cl)
        #:endif
    #:endfor
    #:set _lines = []
    #:set _cur = ''
    #:for _t in _toks
        #:if _cur == ''
            #:set _cur = _t
        #:elif len(_cur) + 1 + len(_t) > width
            #:set _lines = _lines + [_cur + ' &']
            #:set _cur = sentinel + '& ' + _t
        #:else
            #:set _cur = _cur + ' ' + _t
        #:endif
    #:endfor
    #:set _lines = _lines + [_cur]
    $:'\n'.join(_lines)
#:enddef
! New line at end of file is required for FYPP
