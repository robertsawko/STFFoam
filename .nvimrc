"Coding style 
"No tabs, just spaces (Google guide)
set expandtab number ruler
"Indents are four columns (OF guide)

"Line length guard
highlight OverLength ctermbg=red ctermfg=white guibg=#592929
match OverLength /\%81v.\+/

"wildignre allows us to ignore some files when searching and opening from
"within vim
set wildignore+=*.o,*.obj,lnInclude,linux64GccDPOpt,tags,*.dep,gmon.out,
                        \doc,*/build/*,*/run/*,*.swp

let g:ycm_confirm_extra_conf = 0
"Compilation
set makeprg=cd\ build;make

set foldmethod=syntax
set formatprg=clang-format\ 2>/dev/null

"A workaround for slow foldmethod
autocmd InsertEnter * let w:last_fdm=&foldmethod | setlocal foldmethod=manual
autocmd InsertLeave * let &l:foldmethod=w:last_fdm
