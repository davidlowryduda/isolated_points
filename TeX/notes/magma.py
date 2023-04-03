#!/usr/bin/python3.9

import argparse
import sys
import pygments.cmdline as _cmdline

def main(args):
  parser = argparse.ArgumentParser()
  parser.add_argument('-l', dest='lexer', type=str)
  opts, rest = parser.parse_known_args(args[1:])
  if str(opts.lexer).lower() == 'magma':
    args = [__file__, '-l', __file__ + ':Magma', '-x', *rest]
  _cmdline.main(args)

if __name__ == '__main__':
  main(sys.argv)

# ----------------------------------------------------------------------
# -- here's the actual MagmaLexer:

from pygments.lexer import RegexLexer, include, words, bygroups, using, this
from pygments.token import Comment, Operator, Keyword, Name, String, Number,\
                           Punctuation, Whitespace, Text, Generic, Error

class Magma(RegexLexer):
  name = 'Magma'
  aliases = ['magma']
  filenames = ['*.m', '*.mag']

  # the following categories came from an Emacs Magma-mode package
  keywords =       """ _ by default do is select then to where end until
                       catch elif else when case for function if intrinsic
                       procedure repeat try while """
  operators =      """ adj and cat cmpeq cmpne diff div eq ge gt in join
                       le lt meet mod ne notadj notin notsubset or sdiff
                       subset xor not """
  extra_keywords = """ assert assert2 assert3 break clear continue declare delete
                       error error if eval exit
                       forward fprintf freeze iload import load local print
                       printf quit random read readi require requirege
                       requirerange restore return save time vprint vprintf vtime
                       assigned exists forall """
  constructors =   """ car case cop elt ext func hom ideal lideal map
                       ncl pmap proc quo rec recformat rideal sub """
  block_open =     """ try catch case when then else do repeat
                       function procedure intrinsic """
  block_close =    """ until end when elif else catch """

  keywords += operators + extra_keywords + constructors + block_open + block_close
  keywords = keywords.split()

  tokens = {

    'root': [
      include('whitespace-comments'),
      include('error-messages'),
      include('intrinsics'),
      include('number'),
      (r'"', String, 'string'),
      (r'\'[^\']+\'', Name.Builtin),
      (r'\.\.|::|:=|:->|->', Operator),
      (r'[^\w\s]', Punctuation), # do not add +
      (words(keywords, prefix=r'\b', suffix=r'\b'), Keyword.Reserved),
      (r'true|false', Name.Constant),
      (r'[A-Z]+[a-z]+([A-Z]+[a-z]*)?', Name.Builtin),
      (r'\w+', Name),
    ],

    'string': [
      (r'[^\\"]+', String),
      (r'"', String, '#pop'),
      (r'\\.', String.Escape),
    ],

    'whitespace-comments': [
      (r'\s+', Whitespace),
      (r'//.*', Comment),
      (r'/\*', Comment, 'comment-multiline'),
    ],

    'comment-multiline': [
      (r'\*+/', Comment, '#pop'),
      (r'\*+(?!/)', Comment),
      (r'[^*]+', Comment),
    ],

    'number': [
      (r'0[xX][a-fA-F0-9]+', Number.Hex),
      (r'0[bB][01]+', Number.Bin),
      (r'[0-9]+\.([0-9]+|(?!\.))', Number.Float),
      (r'[0-9]+', Number.Integer),
    ],

    'intrinsics': [
      (r'''(?sx)
        (intrinsic) (\s+)
        (\S+?) (\s*)
        ([(].*?[)].*?)
        (-> \s* (?: \S[^,]*? (?: \s* , \s* \S[^,]*?)* \s*))
        (\{)
      ''', bygroups(Keyword, Whitespace, Name.Function, Whitespace,
                    using(this), using(this), Comment), 'comment-intrinsic'),
    ],

    'comment-intrinsic': [
      (r'\{', Comment, '#push'),
      (r'\}', Comment, '#pop'),
      (r'[^\{\}]+', Comment),
    ],

    # This is likely an incomplete list, but who's going to typeset
    # transcripts with errors anyway?
    'error-messages': [
      (r'Runtime error in .*?:.*', Generic.Error),
      (r'Argument types given:', Generic.Error),
      (r'User error:.*', Generic.Error),
      (r'Magma: Internal error', Generic.Error),
      (r'Machine type: \S+', Generic.Error),
      (r'Initial seed: \S+', Generic.Error),
      (r'Time to this point: \S+', Generic.Error),
      (r'Memory usage: \S+?MB', Generic.Error),
      (r'Internal error at \S+, line \S+', Generic.Error),
      (r'Illegal operation', Generic.Error),
      (r'Time:', Generic.Error), # HACK
    ],

  }

# -*- mode: Python; python-indent-offset: 2; python-guess-indent: nil -*-
