/****************************************************************************** 
 * 
 *  file:  ValueArg.h
 * 
 *  Copyright (c) 2003, Michael E. Smoot .
 *  Copyright (c) 2004, Michael E. Smoot, Daniel Aarno.
 *  All rights reverved.
 * 
 *  See the file COPYING in the top directory of this distribution for
 *  more information.
 *  
 *  THE SOFTWARE IS PROVIDED _AS IS_, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 *  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
 *  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 *  DEALINGS IN THE SOFTWARE.  
 *  
 *****************************************************************************/ 


#ifndef TCLAP_SWITCH_ARGUMENT_2_H
#define TCLAP_SWITCH_ARGUMENT_2_H

#include <string>
#include <vector>

#include <tclap/Arg.h>
#include <tclap/Constraint.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define HAVE_SSTREAM
#endif

#if defined(HAVE_SSTREAM)
#include <sstream>
#elif defined(HAVE_STRSTREAM)
#include <strstream>
#else
#error "Need a stringstream (sstream or strstream) to compile!"
#endif

namespace TCLAP {


/**
 * The basic labeled argument that parses a value.
 * This is a template class, which means the type T defines the type
 * that a given object will attempt to parse when the flag/name is matched
 * on the command line.  While there is nothing stopping you from creating
 * an unflagged ValueArg, it is unwise and would cause significant problems.
 * Instead use an UnlabeledValueArg.
 */

class SwitchArg2 : public Arg 
{
    protected:

        /**
         * The value parsed from the command line.
         */
        bool _value;

        /**
         * Extracts the value from the string.
         * Attempts to parse string as type T, if this fails an exception
         * is thrown.
         * \param val - value to be parsed. 
         */
        void _extractValue( const std::string& val );

	public:
 		 /* SwitchArg constructor.
		 * \param flag - The one character flag that identifies this
		 * argument on the command line.
		 * \param name - A one word name for the argument.  Can be
		 * used as a long flag on the command line.
		 * \param desc - A description of what the argument is for or
		 * does.
		 * \param def - The default value for this Switch. 
		 * \param v - An optional visitor.  You probably should not
		 * use this unless you have a very good reason.
		 */
		SwitchArg2(const std::string& flag, 
			      const std::string& name, 
			      const std::string& desc,
			      bool def = false,
				  Visitor* v = NULL);

				  
		/**
		 * SwitchArg constructor.
		 * \param flag - The one character flag that identifies this
		 * argument on the command line.
		 * \param name - A one word name for the argument.  Can be
		 * used as a long flag on the command line.
		 * \param desc - A description of what the argument is for or
		 * does.
		 * \param parser - A CmdLine parser object to add this Arg to
		 * \param def - The default value for this Switch.
		 * \param v - An optional visitor.  You probably should not
		 * use this unless you have a very good reason.
		 */
		SwitchArg2(const std::string& flag, 
			      const std::string& name, 
			      const std::string& desc,
				  CmdLineInterface& parser,
			      bool def = false,
				  Visitor* v = NULL);

        /**
         * Handles the processing of the argument.
         * This re-implements the Arg version of this method to set the
         * _value of the argument appropriately.  It knows the difference
         * between labeled and unlabeled.
         * \param i - Pointer the the current argument in the list.
         * \param args - Mutable list of strings. Passed 
         * in from main().
         */
        virtual bool processArg(int* i, std::vector<std::string>& args); 

        /**
         * Returns the value of the argument.
         */
        bool getValue() ;

        /**
         * Specialization of shortID.
         * \param val - value to be used.
         */
        virtual std::string shortID(const std::string& val = "val") const;

        /**
         * Specialization of longID.
         * \param val - value to be used.
         */
        virtual std::string longID(const std::string& val = "val") const;

};


inline SwitchArg2::SwitchArg2(const std::string& flag, 
	 		         const std::string& name, 
     		   		 const std::string& desc, 
	     	    	 bool _default,
					 Visitor* v )
: Arg(flag, name, desc, false, true, v),
  _value( _default )
{ }

inline SwitchArg2::SwitchArg2(const std::string& flag, 
					const std::string& name, 
					const std::string& desc, 
					CmdLineInterface& parser,
					bool _default,
					Visitor* v )
: Arg(flag, name, desc, false, true, v),
  _value( _default )
{ 
	parser.add( this );
}




/**
 * Implementation of getValue().
 */
inline bool SwitchArg2::getValue() { return _value; }

/**
 * Implementation of processArg().
 */
bool SwitchArg2::processArg(int *i, std::vector<std::string>& args)
{
    if ( _ignoreable && Arg::ignoreRest() )
		return false;

    if ( _hasBlanks( args[*i] ) )
		return false;

    std::string flag = args[*i];

    std::string value = "";
    trimFlag( flag, value );

    if ( argMatches( flag ) )
    {
        if ( _alreadySet )
			throw( CmdLineParseException("Argument already set!", toString()) );

        if ( Arg::delimiter() != ' ' && value == "" )
			throw( ArgParseException( 
							"Couldn't find delimiter for this argument!",
                             toString() ) );

        if ( value == "" )
        {
            (*i)++;
            if ( static_cast<unsigned int>(*i) < args.size() ) 
				_extractValue( args[*i] );
            else
				throw( ArgParseException("Missing a value for this argument!",
                                                    toString() ) );
        }
        else
			_extractValue( value );
				
        _alreadySet = true;
        _checkWithVisitor();
        return true;
    }	
    else
		return false;
}

/**
 * Implementation of shortID.
 */
std::string SwitchArg2::shortID(const std::string& val) const
{
    return Arg::shortID( "0/1" ); 
}

/**
 * Implementation of longID.
 */
std::string SwitchArg2::longID(const std::string& val) const
{
    return Arg::longID( "0/1" );
}

void SwitchArg2::_extractValue( const std::string& val ) 
{
	if ( val=="0" )
		_value=false;
	else if ( val=="1" )
		_value=true;
	else
		throw( ArgParseException("Invalid value '" + val + "'. Value can be only '0' (false) or '1' (true) '" , toString() ) );


}

} // namespace TCLAP

#endif
