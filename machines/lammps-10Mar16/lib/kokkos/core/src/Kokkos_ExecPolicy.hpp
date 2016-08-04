/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXECPOLICY_HPP
#define KOKKOS_EXECPOLICY_HPP

#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_Tags.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Execution policy for work over a range of an integral type.
 *
 * Valid template argument options:
 *
 *  With a specified execution space:
 *    < ExecSpace , WorkTag , { IntConst | IntType } >
 *    < ExecSpace , WorkTag , void >
 *    < ExecSpace , { IntConst | IntType } , void >
 *    < ExecSpace , void , void >
 *
 *  With the default execution space:
 *    < WorkTag , { IntConst | IntType } , void >
 *    < WorkTag , void , void >
 *    < { IntConst | IntType } , void , void >
 *    < void , void , void >
 *
 *  IntType  is a fundamental integral type
 *  IntConst is an Impl::integral_constant< IntType , Blocking >
 *
 *  Blocking is the granularity of partitioning the range among threads.
 */
template< class Arg0 = void , class Arg1 = void , class Arg2 = void 
        , class ExecSpace =
          // The first argument is the execution space,
          // otherwise use the default execution space.
          typename std::conditional
            < Impl::is_execution_space< Arg0 >::value , Arg0
            , Kokkos::DefaultExecutionSpace >::type
        >
class RangePolicy {
private:

  // Default integral type and blocking factor:
  typedef int DefaultIntType ;
  enum { DefaultIntValue = 8 };

  enum { Arg0_Void = Impl::is_same< Arg0 , void >::value };
  enum { Arg1_Void = Impl::is_same< Arg1 , void >::value };
  enum { Arg2_Void = Impl::is_same< Arg2 , void >::value };

  enum { Arg0_ExecSpace = Impl::is_execution_space< Arg0 >::value };

  enum { Arg0_IntConst = Impl::is_integral_constant< Arg0 >::value };
  enum { Arg1_IntConst = Impl::is_integral_constant< Arg1 >::value };
  enum { Arg2_IntConst = Impl::is_integral_constant< Arg2 >::value };

  enum { Arg0_IntType = Impl::is_integral< Arg0 >::value };
  enum { Arg1_IntType = Impl::is_integral< Arg1 >::value };
  enum { Arg2_IntType = Impl::is_integral< Arg2 >::value };

  enum { Arg0_WorkTag = ! Arg0_ExecSpace && ! Arg0_IntConst && ! Arg0_IntType && ! Arg0_Void };
  enum { Arg1_WorkTag =   Arg0_ExecSpace && ! Arg1_IntConst && ! Arg1_IntType && ! Arg1_Void };

  enum { ArgOption_OK = Impl::StaticAssert< (
    ( Arg0_ExecSpace && Arg1_WorkTag && ( Arg2_IntConst || Arg2_IntType ) ) ||
    ( Arg0_ExecSpace && Arg1_WorkTag && Arg2_Void ) ||
    ( Arg0_ExecSpace && ( Arg1_IntConst || Arg1_IntType ) && Arg2_Void ) ||
    ( Arg0_ExecSpace && Arg1_Void && Arg2_Void ) ||
    ( Arg0_WorkTag && ( Arg1_IntConst || Arg1_IntType ) && Arg2_Void ) ||
    ( Arg0_WorkTag && Arg1_Void && Arg2_Void ) ||
    ( ( Arg0_IntConst || Arg0_IntType ) && Arg1_Void && Arg2_Void ) ||
    ( Arg0_Void && Arg1_Void && Arg2_Void )
    ) >::value };

  // The work argument tag is the first or second argument
  typedef typename std::conditional< Arg0_WorkTag , Arg0 ,
          typename std::conditional< Arg1_WorkTag , Arg1 , void
          >::type >::type
    WorkTag ;

  enum { Granularity = Arg0_IntConst ? unsigned(Impl::is_integral_constant<Arg0>::integral_value) : (
                       Arg1_IntConst ? unsigned(Impl::is_integral_constant<Arg1>::integral_value) : (
                       Arg2_IntConst ? unsigned(Impl::is_integral_constant<Arg2>::integral_value) : (
                                       unsigned(DefaultIntValue) ))) };

  // Only accept the integral type if the blocking is a power of two
  static_assert( Impl::is_integral_power_of_two( Granularity )
               , "RangePolicy blocking granularity must be power of two" );

  typedef typename std::conditional< Arg0_IntType , Arg0 ,
          typename std::conditional< Arg1_IntType , Arg1 ,
          typename std::conditional< Arg2_IntType , Arg2 ,
          typename std::conditional< Arg0_IntConst , typename Impl::is_integral_constant<Arg0>::integral_type ,
          typename std::conditional< Arg1_IntConst , typename Impl::is_integral_constant<Arg1>::integral_type ,
          typename std::conditional< Arg2_IntConst , typename Impl::is_integral_constant<Arg2>::integral_type ,
                                                     DefaultIntType
          >::type >::type >::type
          >::type >::type >::type
    IntType ;

  enum { GranularityMask = IntType(Granularity) - 1 };

  ExecSpace m_space ;
  IntType   m_begin ;
  IntType   m_end ;

public:

  //! Tag this class as an execution policy
  typedef ExecSpace    execution_space ;
  typedef RangePolicy  execution_policy ;
  typedef WorkTag      work_tag ;
  typedef IntType      member_type ;

  KOKKOS_INLINE_FUNCTION const execution_space & space() const { return m_space ; }
  KOKKOS_INLINE_FUNCTION member_type begin() const { return m_begin ; }
  KOKKOS_INLINE_FUNCTION member_type end()   const { return m_end ; }

  inline RangePolicy() : m_space(), m_begin(0), m_end(0) {}

  /** \brief  Total range */
  inline
  RangePolicy( const member_type work_begin
             , const member_type work_end
             )
    : m_space()
    , m_begin( work_begin < work_end ? work_begin : 0 )
    , m_end(   work_begin < work_end ? work_end : 0 )
    {}

  /** \brief  Total range */
  inline
  RangePolicy( const execution_space & work_space
             , const member_type work_begin
             , const member_type work_end
             )
    : m_space( work_space )
    , m_begin( work_begin < work_end ? work_begin : 0 )
    , m_end(   work_begin < work_end ? work_end : 0 )
    {}

  /** \brief  Subrange for a partition's rank and size.
   *
   *  Typically used to partition a range over a group of threads.
   */
  struct WorkRange {
    typedef typename RangePolicy::work_tag     work_tag ;
    typedef typename RangePolicy::member_type  member_type ;

    KOKKOS_INLINE_FUNCTION member_type begin() const { return m_begin ; }
    KOKKOS_INLINE_FUNCTION member_type end()   const { return m_end ; }

    /** \brief  Subrange for a partition's rank and size.
     *
     *  Typically used to partition a range over a group of threads.
     */
    KOKKOS_INLINE_FUNCTION
    WorkRange( const RangePolicy & range
             , const int part_rank
             , const int part_size
             )
      : m_begin(0), m_end(0)
      {
        if ( part_size ) {
  
          // Split evenly among partitions, then round up to the granularity.
          const member_type work_part =
            ( ( ( ( range.end() - range.begin() ) + ( part_size - 1 ) ) / part_size )
              + GranularityMask ) & ~member_type(GranularityMask);

          m_begin = range.begin() + work_part * part_rank ;
          m_end   = m_begin       + work_part ;
  
          if ( range.end() < m_begin ) m_begin = range.end() ;
          if ( range.end() < m_end )   m_end   = range.end() ;
        }
      }
  private:
     member_type m_begin ;
     member_type m_end ;
     WorkRange();
     WorkRange & operator = ( const WorkRange & );
  };
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

namespace Experimental {

/** \brief Scratch memory request accepting per team and per thread value
 *
 * An instance of this class can be given as the last argument to a 
 * TeamPolicy constructor. It sets the amount of user requested shared
 * memory for the team.
 */

template< class MemorySpace >
class TeamScratchRequest {
  size_t m_per_team;
  size_t m_per_thread;
  
public:
  TeamScratchRequest(size_t per_team_, size_t per_thread_ = 0):
   m_per_team(per_team_), m_per_thread(per_thread_) {
  } 

  size_t per_team() const {
    return m_per_team;
  }
  size_t per_thread() const {
    return m_per_thread;
  }
  size_t total(const size_t team_size) const {
    return m_per_team + m_per_thread * team_size;
  }
}; 

}

/** \brief  Execution policy for parallel work over a league of teams of threads.
 *
 *  The work functor is called for each thread of each team such that
 *  the team's member threads are guaranteed to be concurrent.
 *
 *  The team's threads have access to team shared scratch memory and
 *  team collective operations.
 *
 *  If the WorkTag is non-void then the first calling argument of the
 *  work functor's parentheses operator is 'const WorkTag &'.
 *  This allows a functor to have multiple work member functions.
 *
 *  template argument option with specified execution space:
 *    < ExecSpace , WorkTag >
 *    < ExecSpace , void >
 *
 *  template argument option with default execution space:
 *    < WorkTag , void >
 *    < void , void >
 */
template< class Arg0 = void
        , class Arg1 = void
        , class ExecSpace =
          // If the first argument is not an execution
          // then use the default execution space.
          typename std::conditional
            < Impl::is_execution_space< Arg0 >::value , Arg0
            , Kokkos::DefaultExecutionSpace >::type
        >
class TeamPolicy {
private:

  enum { Arg0_ExecSpace = Impl::is_execution_space< Arg0 >::value };
  enum { Arg1_Void      = Impl::is_same< Arg1 , void >::value };
  enum { ArgOption_OK   = Impl::StaticAssert< ( Arg0_ExecSpace || Arg1_Void ) >::value };

  typedef typename std::conditional< Arg0_ExecSpace , Arg1 , Arg0 >::type WorkTag ;

public:

  //! Tag this class as an execution policy
  typedef TeamPolicy  execution_policy ;
  typedef ExecSpace   execution_space ;
  typedef WorkTag     work_tag ;

  //----------------------------------------
  /** \brief  Query maximum team size for a given functor.
   *
   *  This size takes into account execution space concurrency limitations and
   *  scratch memory space limitations for reductions, team reduce/scan, and
   *  team shared memory.
   */
  template< class FunctorType >
  static int team_size_max( const FunctorType & );

  /** \brief  Query recommended team size for a given functor.
   *
   *  This size takes into account execution space concurrency limitations and
   *  scratch memory space limitations for reductions, team reduce/scan, and
   *  team shared memory.
   */
  template< class FunctorType >
  static int team_size_recommended( const FunctorType & );

  template< class FunctorType >
  static int team_size_recommended( const FunctorType & , const int&);
  //----------------------------------------
  /** \brief  Construct policy with the given instance of the execution space */
  TeamPolicy( const execution_space & , int league_size_request , int team_size_request , int vector_length_request = 1 );

  TeamPolicy( const execution_space & , int league_size_request , const Kokkos::AUTO_t & , int vector_length_request = 1 );

  /** \brief  Construct policy with the default instance of the execution space */
  TeamPolicy( int league_size_request , int team_size_request , int vector_length_request = 1 );

  TeamPolicy( int league_size_request , const Kokkos::AUTO_t & , int vector_length_request = 1 );

  template<class MemorySpace>
  TeamPolicy( int league_size_request , int team_size_request , const Experimental::TeamScratchRequest<MemorySpace>& team_scratch_memory_request );

  template<class MemorySpace>
  TeamPolicy( int league_size_request , const Kokkos::AUTO_t & , const Experimental::TeamScratchRequest<MemorySpace>& team_scratch_memory_request );

  /** \brief  The actual league size (number of teams) of the policy.
   *
   *  This may be smaller than the requested league size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int league_size() const ;

  /** \brief  The actual team size (number of threads per team) of the policy.
   *
   *  This may be smaller than the requested team size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int team_size() const ;

  /** \brief  Parallel execution of a functor calls the functor once with
   *          each member of the execution policy.
   */
  struct member_type {

    /** \brief  Handle to the currently executing team shared scratch memory */
    KOKKOS_INLINE_FUNCTION
    typename execution_space::scratch_memory_space team_shmem() const ;

    /** \brief  Rank of this team within the league of teams */
    KOKKOS_INLINE_FUNCTION int league_rank() const ;

    /** \brief  Number of teams in the league */
    KOKKOS_INLINE_FUNCTION int league_size() const ;

    /** \brief  Rank of this thread within this team */
    KOKKOS_INLINE_FUNCTION int team_rank() const ;

    /** \brief  Number of threads in this team */
    KOKKOS_INLINE_FUNCTION int team_size() const ;

    /** \brief  Barrier among the threads of this team */
    KOKKOS_INLINE_FUNCTION void team_barrier() const ;

    /** \brief  Intra-team reduction. Returns join of all values of the team members. */
    template< class JoinOp >
    KOKKOS_INLINE_FUNCTION
    typename JoinOp::value_type team_reduce( const typename JoinOp::value_type
                                           , const JoinOp & ) const ;

    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
     *
     *  The highest rank thread can compute the reduction total as
     *    reduction_total = dev.team_scan( value ) + value ;
     */
    template< typename Type >
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const ;

    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
     *          with intra-team non-deterministic ordering accumulation.
     *
     *  The global inter-team accumulation value will, at the end of the
     *  league's parallel execution, be the scan's total.
     *  Parallel execution ordering of the league's teams is non-deterministic.
     *  As such the base value for each team's scan operation is similarly
     *  non-deterministic.
     */
    template< typename Type >
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum ) const ;
  };
};

} // namespace Kokkos

namespace Kokkos {

namespace Impl {

template<typename iType, class TeamMemberType>
struct TeamThreadRangeBoundariesStruct {
private:

  KOKKOS_INLINE_FUNCTION static
  iType ibegin( const iType & arg_begin
              , const iType & arg_end
              , const iType & arg_rank
              , const iType & arg_size
              )
    {
      return arg_begin + ( ( arg_end - arg_begin + arg_size - 1 ) / arg_size ) * arg_rank ;
    }

  KOKKOS_INLINE_FUNCTION static
  iType iend( const iType & arg_begin
            , const iType & arg_end
            , const iType & arg_rank
            , const iType & arg_size
            )
    {
      const iType end_ = arg_begin + ( ( arg_end - arg_begin + arg_size - 1 ) / arg_size ) * ( arg_rank + 1 );
      return end_ < arg_end ? end_ : arg_end ;
    }

public:

  typedef iType index_type;
  const iType start;
  const iType end;
  enum {increment = 1};
  const TeamMemberType& thread;

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct( const TeamMemberType& arg_thread
                                , const iType& arg_end
                                )
    : start( ibegin( 0 , arg_end , arg_thread.team_rank() , arg_thread.team_size() ) )
    , end(   iend(   0 , arg_end , arg_thread.team_rank() , arg_thread.team_size() ) )
    , thread( arg_thread )
    {}

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct( const TeamMemberType& arg_thread
                                , const iType& arg_begin
                                , const iType& arg_end
                                )
    : start( ibegin( arg_begin , arg_end , arg_thread.team_rank() , arg_thread.team_size() ) )
    , end(   iend(   arg_begin , arg_end , arg_thread.team_rank() , arg_thread.team_size() ) )
    , thread( arg_thread )
    {}
};

  template<typename iType, class TeamMemberType>
  struct ThreadVectorRangeBoundariesStruct {
    typedef iType index_type;
    enum {start = 0};
    const iType end;
    enum {increment = 1};

    KOKKOS_INLINE_FUNCTION
    ThreadVectorRangeBoundariesStruct (const TeamMemberType& thread, const iType& count):
      end( count )
    {}
  };

  template<class TeamMemberType>
  struct ThreadSingleStruct {
    const TeamMemberType& team_member;
    KOKKOS_INLINE_FUNCTION
    ThreadSingleStruct(const TeamMemberType& team_member_):team_member(team_member_){}
  };

  template<class TeamMemberType>
  struct VectorSingleStruct {
    const TeamMemberType& team_member;
    KOKKOS_INLINE_FUNCTION
    VectorSingleStruct(const TeamMemberType& team_member_):team_member(team_member_){}
  };
} // namespace Impl

/** \brief  Execution policy for parallel work over a threads within a team.
 *
 *  The range is split over all threads in a team. The Mapping scheme depends on the architecture.
 *  This policy is used together with a parallel pattern as a nested layer within a kernel launched
 *  with the TeamPolicy. This variant expects a single count. So the range is (0,count].
 */
template<typename iType, class TeamMemberType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,TeamMemberType> TeamThreadRange(const TeamMemberType&, const iType& count);

/** \brief  Execution policy for parallel work over a threads within a team.
 *
 *  The range is split over all threads in a team. The Mapping scheme depends on the architecture.
 *  This policy is used together with a parallel pattern as a nested layer within a kernel launched
 *  with the TeamPolicy. This variant expects a begin and end. So the range is (begin,end].
 */
template<typename iType, class TeamMemberType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,TeamMemberType> TeamThreadRange(const TeamMemberType&, const iType& begin, const iType& end);

/** \brief  Execution policy for a vector parallel loop.
 *
 *  The range is split over all vector lanes in a thread. The Mapping scheme depends on the architecture.
 *  This policy is used together with a parallel pattern as a nested layer within a kernel launched
 *  with the TeamPolicy. This variant expects a single count. So the range is (0,count].
 */
template<typename iType, class TeamMemberType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,TeamMemberType> ThreadVectorRange(const TeamMemberType&, const iType& count);

} // namespace Kokkos

#endif /* #define KOKKOS_EXECPOLICY_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

