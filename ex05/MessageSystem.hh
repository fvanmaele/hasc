#ifndef HASC_MESSAGESYSTEM_HH
#define HASC_MESSAGESYSTEM_HH

#include <thread>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <list>
#include <map>
#include <chrono>
#include <iostream>
#include "Barrier.hh"

class MessageSystem {
  int P; // number of threads participating
  Barrier b; // a barrier to sync all threads
  std::mutex mapmx; // for the condition variable
  std::map<std::thread::id,int> myrank; // find rank from thread id
  
  // mailbox part
  struct Envelope {
    int sender; // origin of the message
    int receiver; // destination of the message
    bool read; // set to true by receiver to notify sender
    std::mutex mx; // for the condition variable
    std::condition_variable cv; // for the sender to wait
    Envelope(int _sender) // used for msg list
      : sender(_sender), receiver(-1), read(false)
    {}
    Envelope(int _sender, int _receiver) // used for arecv
      : sender(_sender), receiver(_receiver), read(false)
    {}
    virtual void dummy (void) // enables RTTI
    {}
  };
  template<class T>
  struct TypedEnvelope : public Envelope {
    const T* data;
    T* data_recv;
    TypedEnvelope(int _sender, const T* _data)
      : Envelope(_sender), data(_data)
    {}
    TypedEnvelope(int _sender, int _receiver, T* _data)
      : Envelope(_sender,_receiver), data_recv(_data)
    {}
    virtual void dummy (void)
    {}
  };

  // one message queue for each thread to receive
  using EnvelopeSharedPtr = std::shared_ptr<Envelope>; // internal use
  using MessageID = EnvelopeSharedPtr; // external use
  std::vector<std::list<EnvelopeSharedPtr>> messages;
  // locks to protect the queues
  std::vector<std::mutex> locks;
  // condition variables to wait for a message ariving
  std::vector<std::condition_variable> cvs;

public:
  MessageSystem (int _P) : P(_P), b(_P), locks(_P), messages(_P), cvs(_P)
  {}

  // this must be called once at the start by each thread
  void register_thread (int rank)
  {
    std::unique_lock<std::mutex> ul{mapmx};
    myrank[std::this_thread::get_id()] = rank;
    ul.unlock();
    b.wait(rank);
  }

  // barrier can be called without the rank
  void barrier ()
  {
    b.wait(myrank[std::this_thread::get_id()]);
  }

  // get rank from thread id
  int get_rank ()
  {
    return myrank[std::this_thread::get_id()];
  }

  // blocking typesafe send
  template<typename T>
  void send (int to_rank, const T& data)
  {
    int me = myrank[std::this_thread::get_id()];
    
    // check if message goes to itself
    if (me==to_rank)
      {
	std::cout << "thread " << to_rank << " sends message to itself" << std::endl;
	return;
      }
    
    // make envelope
    using TEnvelope = TypedEnvelope<T>;
    EnvelopeSharedPtr penv(new TEnvelope(me,&data));
    
    // put envelope in message list of receiving thread
    locks[to_rank].lock();
    messages[to_rank].push_back(penv); // append element at the end
    locks[to_rank].unlock();
        
    // notify receiving thread that new message has arrived
    cvs[to_rank].notify_one();

    // block on the envelopes flag until message has been read
    std::unique_lock<std::mutex> ul{penv->mx};
    penv->cv.wait(ul,[&]{return penv->read;});
    // shared_ptr to envelope has been removed from list by receiver
    // so we can destroy the envelope now
  }

  // blocking receive
  template<typename T>
  void recv (int from_rank, T& data)
  {
    int me = myrank[std::this_thread::get_id()];

    // check if try to receive message from myself
    if (me==from_rank)
      {
	std::cout << "thread " << from_rank << " tries to receive message from itself" << std::endl;
	return;
      }

    // wait for a message arriving from the given rank
    while (true) // hopefully we exit again ...
      {
	// grab lock protecting message list
	std::unique_lock<std::mutex> ul{locks[me]};

	// wait for the message to arrive
	cvs[me].wait(ul,[&]{
	    for (auto iter=messages[me].begin(); iter!=messages[me].end(); ++iter)
	      if ((**iter).sender==from_rank) return true;
	    return false;});

	// we have the lock and the message;
	// find first entry in list from sender
	auto it = messages[me].end();
	for (auto iter=messages[me].begin(); iter!=messages[me].end(); ++iter)
	  if ((**iter).sender==from_rank)
	    {
	      it = iter;
	      break;
	    }
	// now really receive the message; still in the lock
	if (it!=messages[me].end())
	  {
	    // type safe downcast of Envelope* to TypedEnvelope<T>*
	    TypedEnvelope<T>* p = dynamic_cast<TypedEnvelope<T>*>(it->get()); 
	    if (p!=0)
	      data = *(p->data); // now thats the copy using an assignement operator
	    else
	      std::cout << "type mismatch in recv" << std::endl;
	    // signal to sender that message hes been received
	    std::unique_lock<std::mutex> ulp{(**it).mx};
	    (**it).read = true;
	    (**it).cv.notify_one();
	    ulp.unlock();
	    // erase the element from the list; we still haave the lock ul!
	    messages[me].erase(it);
	    // and leave while loop
	    break; // releases the lock ul as well!
	  }
      }
  }

  // blocking receive from any
  // reports from whom the message was received
  template<typename T>
  void recvany (int& from_rank, T& data)
  {
    int me = myrank[std::this_thread::get_id()];

    // wait for a message to arrive
    while (true) // hopefully we exit again ...
      {
	// grab lock protecting message list
	std::unique_lock<std::mutex> ul{locks[me]};

	// wait for a message to arrive
	cvs[me].wait(ul,[&]{return messages[me].size()>0;});

	// we have the lock and the message;
	// get first entry in the list
	auto it = messages[me].begin();

	// now really receive the message; still in the lock
	if (it!=messages[me].end())
	  {
	    // type safe downcast of Envelope* to TypedEnvelope<T>*
	    TypedEnvelope<T>* p = dynamic_cast<TypedEnvelope<T>*>(it->get()); 
	    if (p!=0)
	      {
		from_rank = p->sender; // get the sender
		data = *(p->data); // now thats the copy using an assignement operator;
	      }
	    else
	      std::cout << "type mismatch in recv" << std::endl;
	    // signal to sender that message has been received
	    std::unique_lock<std::mutex> ulp{(**it).mx};
	    (**it).read = true;
	    (**it).cv.notify_one();
	    ulp.unlock();
	    // erase the element from the list
	    messages[me].erase(it);
	    // and leave while loop
	    break; // releases the lock as well!
	  }
      }
  }

  // returns true if a message from the given rank is available
  bool rprobe (int from_rank)
  {
    int me = myrank[std::this_thread::get_id()];
    
    // check if try to receive message from myself
    if (me==from_rank)
      {
	std::cout << "thread " << from_rank << " tries to probe from itself" << std::endl;
	return false;
      }

    // grab lock protecting message list
    std::unique_lock<std::mutex> ul{locks[me]};

    // scan the message list
    for (auto iter=messages[me].begin(); iter!=messages[me].end(); ++iter)
      if ((**iter).sender==from_rank)
	return true;
    return false;
  }

  // nonblocking typesafe send
  template<typename T>
  MessageID asend (int to_rank, const T& data)
  {
    int me = myrank[std::this_thread::get_id()];
    
    // check if message goes to itself
    if (me==to_rank)
      {
	std::cout << "thread " << to_rank << " sends message to itself" << std::endl;
      }
    
    // make envelope which is also the message id
    using TEnvelope = TypedEnvelope<T>;
    EnvelopeSharedPtr penv(new TEnvelope(me,&data));
    
    // put envelope in message list of receiving thread
    locks[to_rank].lock();
    messages[to_rank].push_back(penv); // append element at the end
    locks[to_rank].unlock();
        
    // notify receiving thread that new message has arrived
    cvs[to_rank].notify_one();

    return penv;
  }

  // test whether message has been transfered
  template<typename T>
  bool test (MessageID& msg)
  {
    if (msg->receiver==-1)
      {
	// test if asend succeeded
	std::lock_guard<std::mutex> lg{msg->mx};
	return msg->read;
      }
    else
      {
	// test if arecv succeeded
	int me = myrank[std::this_thread::get_id()];

	// grab lock protecting message list
	std::unique_lock<std::mutex> ul{locks[me]};
	
	// sscan the list to see if message arrived
	auto it = messages[me].end();
	for (auto iter=messages[me].begin(); iter!=messages[me].end(); ++iter)
	  if ((**iter).sender==msg->sender)
	    {
	      it = iter;
	      break;
	    }
	// return if message did not arrive yet
	if (it==messages[me].end()) return false;

	// now really receive the message; still in the lock
	// type safe downcast of Envelope* to TypedEnvelope<T>*
	TypedEnvelope<T>* psend = dynamic_cast<TypedEnvelope<T>*>(it->get());
	TypedEnvelope<T>* precv = dynamic_cast<TypedEnvelope<T>*>(msg.get());
	if (psend!=0 && precv!=0)
	  *(precv->data_recv) = *(psend->data); // now thats the copy using an assignement operator
	else
	  std::cout << "type mismatch in recv" << std::endl;
	// signal to sender that message hes been received
	// this would even work with a blocking send!
	std::unique_lock<std::mutex> ulp{(**it).mx};
	(**it).read = true;
	(**it).cv.notify_one();
	ulp.unlock();
	// erase the element from the list; we still haave the lock ul!
	messages[me].erase(it);
	return true;
      }
  }

  // nonblocking receive
  template<typename T>
  MessageID arecv (int from_rank, T& data)
  {
    int to_me = myrank[std::this_thread::get_id()];

    // check if try to receive message from myself
    if (to_me==from_rank)
      {
	std::cout << "thread " << from_rank << " tries to receive message from itself" << std::endl;
      }

    // make receiver envelope which is also the message id
    using TEnvelope = TypedEnvelope<T>;
    EnvelopeSharedPtr penv(new TEnvelope(from_rank,to_me,&data));

    return penv;
  }
};

// here is the global data
std::shared_ptr<MessageSystem> ms(nullptr);

// call this function once at the start before any threads
void initialize_message_system (int P)
{
  std::atomic_thread_fence(std::memory_order_seq_cst);
  ms = std::make_shared<MessageSystem>(P);
  std::atomic_thread_fence(std::memory_order_seq_cst);
}

// this function must be called once by each thread
void register_thread (int rank)
{
  ms->register_thread(rank);
}

// and free standing functions for send and recv

// blocking send
template<typename T>
void send (int to_rank, const T& data)
{
  ms->send(to_rank,data);
}

// blocking receive
template<typename T>
void recv (int from_rank, T& data)
{
  ms->recv(from_rank,data);
}

// blocking receive from any
template<typename T>
void recvany (int& from_rank, T& data)
{
  ms->recvany(from_rank,data);
}

// nonblocking probe if receive would block
// returns true if a message from the given rank is available
bool rprobe (int from_rank)
{
  return ms->rprobe(from_rank);
}

// wait for all
void barrier ()
{
  return ms->barrier();
}

// nonblocking send
template<typename T>
MessageSystem::MessageID asend (int to_rank, const T& data)
{
  return ms->asend(to_rank,data);
}

// nonblocking receive
template<typename T>
MessageSystem::MessageID arecv (int from_rank, T& data)
{
  return ms->arecv(from_rank,data);
}

// nonblocking test
template<typename T>
bool test (MessageSystem::MessageID& msg)
{
  return ms->test<T>(msg);
}

#endif
